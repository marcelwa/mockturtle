/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file rank_view.hpp
  \brief Implements rank orders for a network

  \author Marcel Walter
*/

#pragma once

#include "../traits.hpp"
#include "../utils/node_map.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

namespace mockturtle
{

/*! \brief Implements rank orders for a network.
 *
 * This view assigns manipulable relative orders to the nodes of each level.
 * A sequence of nodes in the same level is called a rank. The width
 * of a rank is thereby defined as the number of nodes in the rank. The width of
 * the network is equal to the width of the widest rank. This view implements
 * functions to retrieve the assigned node position within a rank (`rank_position`)
 * as well as to fetch the node in a certain rank at a certain position (`at_rank_position`).
 * The ranks are assigned at construction and can be manipulated by calling the `swap` function.
 *
 * This view also automatically inserts new nodes into their respective rank (at the end),
 * however, it does not update the information, when modifying or deleting nodes.
 *
 * **Required network functions:**
 * -  foreach_node
 * -  is_constant
 * -  level
 * -  depth
 *
 * Example
 *
   \verbatim embed:rst

   .. code-block:: c++

      // create network somehow
      aig_network aig = ...;

      // create a depth view on the network
      depth_view aig_depth{aig};
      // create a rank view on the depth view
      rank_view aig_rank{aig_depth};

      // print width
      std::cout << "Width: " << aig_rank.width() << "\n";
   \endverbatim
 */
template<class Ntk, bool has_rank_interface = has_rank_position_v<Ntk>&& has_at_rank_position_v<Ntk>&& has_width_v<Ntk>>
class rank_view
{
};

template<class Ntk>
class rank_view<Ntk, true> : public Ntk
{
public:
  rank_view( Ntk const& ntk ) : Ntk( ntk )
  {
  }
};

template<class Ntk>
class rank_view<Ntk, false> : public Ntk
{
public:
  using storage = typename Ntk::storage;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  explicit rank_view()
      : Ntk(), rank_pos{ *this }, ranks{}, max_rank_width{ 0 }
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
    static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
    static_assert( has_level_v<Ntk>, "Ntk does not implement the level method" );
    static_assert( has_depth_v<Ntk>, "Ntk does not implement the depth method" );

    add_event = Ntk::events().register_add_event( [this]( auto const& n )
                                                  { on_add( n ); } );
  }

  /*! \brief Standard constructor.
   *
   * \param ntk Base network
   */
  explicit rank_view( Ntk const& ntk )
      : Ntk{ ntk }, rank_pos{ ntk }, ranks{ ntk.depth() + 1 }, max_rank_width{ 0 }
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
    static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
    static_assert( has_level_v<Ntk>, "Ntk does not implement the level method" );
    static_assert( has_depth_v<Ntk>, "Ntk does not implement the depth method" );

    init_ranks();

    add_event = Ntk::events().register_add_event( [this]( auto const& n )
                                                  { on_add( n ); } );
  }

  /*! \brief Copy constructor. */
  explicit rank_view( rank_view<Ntk, false> const& other )
      : Ntk{ other }, rank_pos{ other.rank_pos }, ranks{ other.ranks }, max_rank_width{ other.max_rank_width }
  {
    add_event = Ntk::events().register_add_event( [this]( auto const& n )
                                                  { on_add( n ); } );
  }

  rank_view<Ntk, false>& operator=( rank_view<Ntk, false> const& other )
  {
    /* delete the event of this network */
    Ntk::events().release_add_event( add_event );

    /* update the base class */
    this->_storage = other._storage;
    this->_events = other._events;

    /* copy */
    rank_pos = other.rank_pos;
    ranks = other.ranks;
    max_rank_width = other.max_rank_width;

    /* register new event in the other network */
    add_event = Ntk::events().register_add_event( [this]( auto const& n )
                                                  { on_add( n ); } );

    return *this;
  }

  ~rank_view()
  {
    Ntk::events().release_add_event( add_event );
  }

  uint32_t rank_position( node const& n ) const noexcept
  {
    assert( !this->is_constant( n ) && "node must not be constant" );

    return rank_pos[n];
  }

  node at_rank_position( uint32_t const level, uint32_t const pos ) const noexcept
  {
    assert( level < ranks.size() && "level must be less than the number of ranks" );
    assert( pos < ranks[level].size() && "pos must be less than the number of nodes in rank" );

    return ranks[level][pos];
  }

  uint32_t width() const noexcept
  {
    return max_rank_width;
  }

  void swap( node const& n1, node const& n2 ) noexcept
  {
    assert( this->level( n1 ) == this->level( n2 ) && "nodes must be in the same rank" );

    auto &pos1 = rank_pos[n1], pos2 = rank_pos[n2];

    std::swap( ranks[this->level( n1 )][pos1], ranks[this->level( n2 )][pos2] );
    std::swap( pos1, pos2 );
  }
  /**
   * Overrides the base class method to also call the add_event on create_pi().
   *
   * @note This can (and in fact will) lead to issues if Ntk already calls add_event functions on create_pi()!
   *
   * @return Newly created PI signal.
   */
  signal create_pi()
  {
    auto const n = Ntk::create_pi();
    this->resize_levels(); // this line assumes a depth_view for Ntk
    on_add( this->get_node( n ) );
    return n;
  }

private:
  node_map<uint32_t, Ntk> rank_pos;
  std::vector<std::vector<node>> ranks;
  uint32_t max_rank_width;

  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;

  void insert_in_rank( node const& n ) noexcept
  {
    auto& rank = ranks[this->level( n )];
    rank_pos[n] = rank.size();
    rank.push_back( n );
    max_rank_width = std::max( max_rank_width, static_cast<uint32_t>( rank.size() ) );
  }

  void on_add( node const& n ) noexcept
  {
    if ( this->level( n ) >= ranks.size() )
    {
      // add sufficient ranks to store the new node
      ranks.insert( ranks.end(), this->level( n ) - ranks.size() + 1, {} );
    }
    rank_pos.resize();

    insert_in_rank( n );
  }

  void init_ranks() noexcept
  {
    this->foreach_node( [this]( auto const& n )
                        {
                          if (!this->is_constant(n))
                          {
                            insert_in_rank(n);
                          } } );
  }
};

template<class T>
rank_view( T const& ) -> rank_view<T>;

} // namespace mockturtle
