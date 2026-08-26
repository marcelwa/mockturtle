/*!
  \file acd.hpp
  \brief Ashenhurst-Curtis decomposition as a resynthesis function

  This is the non-disjoint counterpart of `dsd_resynthesis` (dsd.hpp), and it is written
  deliberately in that file's shape -- decompose, then hand each irreducible block to a
  nested resynthesis function -- so that the two can be A/B-ed as drop-in replacements for
  one another inside `node_resynthesis`, `cut_rewriting` and `refactoring`.

  Where `dsd_decomposition` peels off single variables and stops at a *prime* remainder,
  Ashenhurst-Curtis decomposition splits the support into a free set and a bound set that
  may **share** variables, and emits a cascade of `lut_size`-input blocks.  Every block is
  then realised in `Ntk` by the nested resynthesis function.  The decomposer is ABC's
  (`Tempia Calvino`), vendored under `algorithms/acd/`; `acd_iface::decompose` returns the
  cascade in ABC's `decompArray` byte layout and the decoder in `acd_expand.hpp` is reused
  verbatim to read it.

  Scope, stated because it is easy to get wrong:

  - The two-level engine accepts at most **11 variables** (`ac_decomposition.hpp`,
    `max_num_vars`), and in practice accepts ~100 % of functions up to 10 variables, ~69 %
    at 11 and nothing above -- measured, see `2026-08-26-acd-resyn-census`.  Outside that
    band the functor emits nothing, unless `use_fallback` is set (read its documentation
    before setting it).
  - The blocks the cascade emits have up to `lut_size` inputs, so the nested function must
    cover that width.  `lut_size = 4` with `xag_npn_resynthesis` is exact and fast;
    `lut_size = 6` needs something that covers six variables (`exact_resynthesis`,
    `sop_factoring`, or a nested `dsd_resynthesis`).
  - ACD declines a lot at the top of its range, so wrapping this in
    `cached_resynthesis` (cached.hpp) is recommended: its *blacklist* half is the one that
    pays here.

  \author agentic-synthesis
*/

#pragma once

#include <cstdint>
#include <vector>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../../traits.hpp"
#include "../acd/acd_wrapper.hpp"
#include "../acd_expand.hpp"
#include "traits.hpp"

namespace mockturtle
{

/*! \brief Parameters for `acd_resynthesis`. */
struct acd_resynthesis_params
{
  /*! \brief Input size of the blocks the cascade is built from.
   *
   * Must be at least 2 and at most 6.  The nested resynthesis function has to be able to
   * realise an arbitrary function of this many variables.
   */
  uint32_t lut_size{ 4u };

  /*! \brief Do not call the decomposer on functions with more than this many variables.
   *
   * The engine's own hard ceiling is 11; there is no point asking above it and every
   * question costs a truth-table copy.  Lowering it to 10 avoids the band in which ACD
   * refuses roughly a third of what it is offered.
   */
  uint32_t max_num_vars{ 11u };

  /*! \brief Hand the whole function to the nested engine when ACD declines.
   *
   * Off by default, and that is deliberate.  The nested engine's job is to realise the
   * cascade's *blocks*, which have at most `lut_size` inputs; a decline hands it the whole
   * function instead, which may have up to 11 variables.  `xag_npn_resynthesis` reads past
   * its 4-input database if asked that (a segfault, not an assertion), so enabling this is
   * only safe with an engine of unbounded width such as `sop_factoring`.
   *
   * With it off, a decline emits nothing.  `cut_rewriting` and `refactoring` treat that as
   * "no candidate for this cut" and keep the original structure, which is the behaviour a
   * *candidate-offering* operator wants.  `node_resynthesis` needs total coverage, so
   * compose this with a second engine there (see `composed.hpp`) rather than setting this.
   */
  bool use_fallback{ false };
};

/*! \brief Statistics of `acd_resynthesis`. */
struct acd_resynthesis_stats
{
  /*! \brief Calls in which ACD produced a cascade. */
  uint32_t num_accepted{ 0u };
  /*! \brief Calls ACD declined (out of range, or no decomposition). */
  uint32_t num_declined{ 0u };
  /*! \brief Calls delegated to the fallback because the function was narrow enough. */
  uint32_t num_trivial{ 0u };
  /*! \brief Blocks emitted in total. */
  uint32_t num_blocks{ 0u };
  /*! \brief Calls in which the nested function failed to realise a block. */
  uint32_t num_block_failures{ 0u };
  /*! \brief Widest block the cascade ever asked the nested function for. */
  uint32_t max_block_fanins{ 0u };
  /*! \brief Fanin count of the first block the nested function could not realise. */
  uint32_t first_failed_fanins{ 0u };
};

/*! \brief Resynthesis function based on Ashenhurst-Curtis decomposition.
 *
 * Can be passed to `node_resynthesis`, `cut_rewriting` and `refactoring`.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      xag_network xag = ...;

      xag_npn_resynthesis<xag_network> blocks;            // realises the <= 4-input blocks
      acd_resynthesis<xag_network, decltype( blocks )> resyn( blocks );
      cut_rewriting( xag, resyn );
      xag = cleanup_dangling( xag );
   \endverbatim
 */
template<class Ntk, class ResynthesisFn>
class acd_resynthesis
{
public:
  explicit acd_resynthesis( ResynthesisFn& resyn_fn, acd_resynthesis_params const& ps = {},
                            acd_resynthesis_stats* pst = nullptr )
      : _resyn_fn( resyn_fn ), _ps( ps ), _pst( pst )
  {
    assert( _ps.lut_size >= 2u && _ps.lut_size <= 6u );
  }

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function,
                   LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    std::vector<signal<Ntk>> leaves( begin, end );
    uint32_t const num_vars = function.num_vars();

    /* narrow enough to be one block: this is the fallback's job, not ours */
    if ( num_vars <= _ps.lut_size )
    {
      bump( &acd_resynthesis_stats::num_trivial );
      if ( _ps.use_fallback )
        call_fallback( ntk, function, leaves, fn );
      return;
    }

    if ( num_vars > _ps.max_num_vars || num_vars > 11u )
    {
      bump( &acd_resynthesis_stats::num_declined );
      if ( _ps.use_fallback )
        call_fallback( ntk, function, leaves, fn );
      return;
    }

    /* `decompose` mutates its input, so hand it a copy */
    std::vector<uint64_t> bits( function._bits.begin(), function._bits.end() );
    unsigned char arr[256];
    if ( acd_iface::decompose( bits.data(), num_vars, _ps.lut_size, arr ) != 0 )
    {
      bump( &acd_resynthesis_stats::num_declined );
      if ( _ps.use_fallback )
        call_fallback( ntk, function, leaves, fn );
      return;
    }

    /* ---- walk ABC's decompArray and rebuild the cascade in `Ntk` --------------
     * layout: [0] byte count, [1] number of blocks (the last one is the root); then per
     * block: fanin count, that many fanin indices (< num_vars: a leaf of this call;
     * otherwise index - num_vars into the blocks already built), then the block's truth
     * table little-endian, one byte at a time.  Decoded by `detail::acd_decode_tt`. */
    uint32_t const num_blocks = arr[1];
    std::vector<signal<Ntk>> blocks;
    blocks.reserve( num_blocks );

    uint32_t byte_p = 2u;
    bool ok = true;
    for ( uint32_t i = 0u; i < num_blocks && ok; ++i )
    {
      uint32_t const num_fanins = arr[byte_p++];
      std::vector<signal<Ntk>> block_leaves;
      block_leaves.reserve( num_fanins );
      for ( uint32_t j = 0u; j < num_fanins; ++j )
      {
        uint32_t const idx = arr[byte_p++];
        block_leaves.push_back( idx < num_vars ? leaves[idx] : blocks[idx - num_vars] );
      }
      if ( _pst != nullptr && num_fanins > _pst->max_block_fanins )
        _pst->max_block_fanins = num_fanins;
      auto const btt = detail::acd_decode_tt( arr, byte_p, num_fanins );

      bool got = false;
      signal<Ntk> s = ntk.get_constant( false );
      auto const on_signal = [&]( signal<Ntk> const& candidate ) {
        if ( !got )
        {
          s = candidate;
          got = true;
        }
        return true; /* keep the first; the caller ranks whole candidates, not blocks */
      };

      if constexpr ( has_set_bounds_v<ResynthesisFn> )
      {
        _resyn_fn.set_bounds( num_fanins, std::nullopt );
      }
      _resyn_fn( ntk, btt, block_leaves.begin(), block_leaves.end(), on_signal );

      if ( !got )
      {
        /* the nested engine could not realise a block: abandon this candidate entirely
         * rather than emit a partial structure.  Nodes already created are dangling and
         * are removed by the usual `cleanup_dangling`. */
        ok = false;
        bump( &acd_resynthesis_stats::num_block_failures );
        if ( _pst != nullptr && _pst->first_failed_fanins == 0u )
          _pst->first_failed_fanins = num_fanins;
        break;
      }
      blocks.push_back( s );
    }

    if ( !ok )
    {
      bump( &acd_resynthesis_stats::num_declined );
      if ( _ps.use_fallback )
        call_fallback( ntk, function, leaves, fn );
      return;
    }

    if ( _pst != nullptr )
    {
      ++_pst->num_accepted;
      _pst->num_blocks += num_blocks;
    }
    fn( blocks.back() );
  }

  void clear_functions()
  {
    if constexpr ( has_clear_functions_v<ResynthesisFn> )
    {
      _resyn_fn.clear_functions();
    }
  }

  void add_function( signal<Ntk> const& s, kitty::dynamic_truth_table const& tt )
  {
    if constexpr ( has_add_function_v<ResynthesisFn, Ntk> )
    {
      _resyn_fn.add_function( s, tt );
    }
  }

private:
  template<typename Fn>
  void call_fallback( Ntk& ntk, kitty::dynamic_truth_table const& function,
                      std::vector<signal<Ntk>>& leaves, Fn&& fn ) const
  {
    if constexpr ( has_set_bounds_v<ResynthesisFn> )
    {
      _resyn_fn.set_bounds( static_cast<uint32_t>( leaves.size() ), std::nullopt );
    }
    _resyn_fn( ntk, function, leaves.begin(), leaves.end(), fn );
  }

  void bump( uint32_t acd_resynthesis_stats::*field ) const
  {
    if ( _pst != nullptr )
      ++( _pst->*field );
  }

  ResynthesisFn& _resyn_fn;
  acd_resynthesis_params _ps;
  acd_resynthesis_stats* _pst;
};

} /* namespace mockturtle */
