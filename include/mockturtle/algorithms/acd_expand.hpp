/*!
  \file acd_expand.hpp
  \brief Materialises wide LUTs as Ashenhurst-Curtis decomposed LUT structures

  A k-LUT network produced by `lut_map` with `acd_lut_size != 0` may contain nodes with
  more than `acd_lut_size` fanins.  Those nodes were costed during mapping by the number
  of `acd_lut_size`-LUTs an Ashenhurst-Curtis decomposition of their function requires;
  this pass performs that decomposition and replaces each such node by the corresponding
  structure, so that the resulting network is a legal `acd_lut_size`-LUT netlist.

  It is the mockturtle counterpart of ABC's `Abc_DecRecordToHop`.
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../networks/klut.hpp"
#include "../traits.hpp"
#include "../views/topo_view.hpp"
#include "acd/acd_wrapper.hpp"

namespace mockturtle
{

namespace detail
{

/*! \brief Decodes one LUT truth table from ABC's decomposition-array byte layout. */
inline kitty::dynamic_truth_table acd_decode_tt( unsigned char const* arr, uint32_t& byte_p, uint32_t num_fanins )
{
  uint32_t const num_words = ( num_fanins <= 6 ) ? 1u : ( 1u << ( num_fanins - 6 ) );
  uint32_t const num_bytes = ( num_fanins <= 3 ) ? 1u : ( 1u << ( std::min<uint32_t>( num_fanins, 6 ) - 3 ) );

  std::vector<uint64_t> words( num_words, 0u );
  for ( uint32_t j = 0; j < num_words; ++j )
  {
    for ( uint32_t k = 0; k < num_bytes; ++k )
      words[j] |= static_cast<uint64_t>( arr[byte_p++] ) << ( k << 3 );
  }

  /* replicate the sub-word pattern the way ABC's reader does */
  if ( num_fanins == 2 )
    words[0] |= words[0] << 4;
  uint32_t nb = num_bytes;
  while ( nb < 4 )
  {
    words[0] |= words[0] << ( nb << 3 );
    nb <<= 1;
  }

  kitty::dynamic_truth_table tt( num_fanins );
  for ( uint32_t j = 0; j < num_words; ++j )
    tt._bits[j] = words[j];
  tt.mask_bits();
  return tt;
}

} /* namespace detail */

/*! \brief Statistics of `acd_expand`. */
struct acd_expand_stats
{
  /*! \brief Number of wide LUTs that were decomposed. */
  uint32_t num_decomposed{ 0u };
  /*! \brief Number of LUTs the decomposed structures contain in total. */
  uint32_t num_luts_created{ 0u };
  /*! \brief Number of wide LUTs for which no decomposition was found (a bug if non-zero). */
  uint32_t num_failed{ 0u };
};

/*! \brief Replaces every LUT wider than `lut_size` by its AC decomposition.
 *
 * \param ntk k-LUT network, possibly containing nodes with up to 11 fanins
 * \param lut_size target LUT size
 * \return an equivalent k-LUT network in which no node has more than `lut_size` fanins
 */
inline klut_network acd_expand( klut_network const& ntk, uint32_t lut_size, acd_expand_stats* pst = nullptr )
{
  acd_expand_stats st;
  klut_network res;

  std::vector<klut_network::signal> old_to_new( ntk.size() );
  old_to_new[ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) )] = res.get_constant( false );
  if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    old_to_new[ntk.node_to_index( ntk.get_node( ntk.get_constant( true ) ) )] = res.get_constant( true );

  ntk.foreach_pi( [&]( auto const& n ) {
    old_to_new[ntk.node_to_index( n )] = res.create_pi();
  } );

  topo_view<klut_network> topo{ ntk };
  topo.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return;

    std::vector<klut_network::signal> children;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      children.push_back( old_to_new[ntk.node_to_index( f )] );
    } );

    auto tt = ntk.node_function( n );

    if ( children.size() <= lut_size )
    {
      old_to_new[ntk.node_to_index( n )] = res.create_node( children, tt );
      return;
    }

    /* wide node: decompose */
    unsigned char arr[256];
    std::vector<uint64_t> bits( tt._bits.begin(), tt._bits.end() );
    if ( acd_iface::decompose( bits.data(), static_cast<uint32_t>( children.size() ), lut_size, arr ) != 0 )
    {
      /* should not happen: the mapper only kept this cut because ACD said it was feasible */
      ++st.num_failed;
      old_to_new[ntk.node_to_index( n )] = res.create_node( children, tt );
      return;
    }

    ++st.num_decomposed;
    uint32_t const num_luts = arr[1];
    st.num_luts_created += num_luts;

    std::vector<klut_network::signal> structure; /* signals of the intermediate LUTs */
    uint32_t byte_p = 2;
    klut_network::signal root = res.get_constant( false );
    for ( uint32_t i = 0; i < num_luts; ++i )
    {
      uint32_t const num_fanins = arr[byte_p++];
      std::vector<klut_network::signal> lut_children;
      for ( uint32_t j = 0; j < num_fanins; ++j )
      {
        uint32_t const idx = arr[byte_p++];
        if ( idx < children.size() )
          lut_children.push_back( children[idx] );
        else
          lut_children.push_back( structure[idx - children.size()] );
      }
      auto const ltt = detail::acd_decode_tt( arr, byte_p, num_fanins );
      root = res.create_node( lut_children, ltt );
      structure.push_back( root );
    }

    old_to_new[ntk.node_to_index( n )] = root;
  } );

  ntk.foreach_po( [&]( auto const& f ) {
    res.create_po( old_to_new[ntk.node_to_index( f )] );
  } );

  if ( pst != nullptr )
    *pst = st;

  return res;
}

} /* namespace mockturtle */
