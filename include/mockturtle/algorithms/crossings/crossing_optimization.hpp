/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022 EPFL
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
\file crossing_optimization.hpp
\brief reorder the nodes of a network to reduce the number of crossings
\author Marcel Walter
*/

#pragma once

#include "../../traits.hpp"
#include "../../utils/hash_functions.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../views/topo_view.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace mockturtle
{

struct crossing_optimization_params
{
  // number of sweeps
  uint16_t rho{ 5 };
};

struct crossing_optimization_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };
  /*! \brief Number of crossings */
  uint32_t num_crossings{ 0 };
};

namespace detail
{

template<typename Ntk>
class crossing_optimization_impl
{
public:
  crossing_optimization_impl( Ntk& src, crossing_optimization_params const& ps, crossing_optimization_stats& st ) noexcept
      : ntk{ src },
        params{ ps },
        stats{ st },
        get_block{ ntk }
  {
  }

  void run() noexcept
  {
    stopwatch t{ stats.time_total };

    initialize_block_list();

    global_sifting();

    // sort each rank by the block list order
    for ( auto level = 0u; level < ntk.depth(); ++level )
    {
      ntk.sort_rank( level, [this]( const auto& lhs, const auto& rhs )
                     { return get_block[lhs]->pi < get_block[rhs]->pi; } );
    }
  }

private:
  /**
   * Network.
   */
  Ntk& ntk;
  /**
   * Parameters.
   */
  crossing_optimization_params const& params;
  /**
   * Statistics.
   */
  crossing_optimization_stats& stats;

  struct block
  {
    // nodes belonging to this block
    std::vector<node<Ntk>> nodes{};

    // fanouts and fanins
    std::vector<node<Ntk>> N_plus{}, N_minus{};
    // index pointers to fanouts and fanins
    std::vector<std::size_t> I_plus{}, I_minus{};

    // position in global ordering
    std::size_t pi{};

    // returns the node with the lowest level in the block
    node<Ntk> upper() const noexcept
    {
      return nodes.front();
    }

    // returns the node with the highest level in the block
    node<Ntk> lower() const noexcept
    {
      return nodes.back();
    }

    std::pair<uint32_t, uint32_t> span;
  };

  // assigns a block to each node
  node_map<std::shared_ptr<block>, Ntk> get_block;

  using block_list = std::vector<std::shared_ptr<block>>;
  block_list B{};

  // gathers nodes that could belong to the same block (1-in-1-out chains)
  template<typename TopoNtk>
  void gather_block_nodes( const TopoNtk& topo_ntk, const node<TopoNtk>& n, std::vector<node<TopoNtk>>& nodes ) noexcept
  {
    // if node n was already visited, skip it
    if ( topo_ntk.visited( n ) == topo_ntk.trav_id() )
    {
      return;
    }

    // otherwise, mark it as visited and add it to the node list
    topo_ntk.set_visited( n, topo_ntk.trav_id() );
    nodes.push_back( n );

    // if node n has only one input and one output, it can be part of a "dummy node block"
    if ( topo_ntk.fanin_size( n ) == 1 && topo_ntk.fanout_size( n ) == 1 )
    {
      topo_ntk.foreach_fanout( n, [this, &topo_ntk, &nodes]( const auto& fon )
                               { gather_block_nodes( topo_ntk, fon, nodes ); } );
    }
  }

  // computes the level span of a block given by the upper and lower levels it encompasses; if the block is just a single node (with no dummy vertices), the level span is 0
  // should a node's fanouts skip any levels, the lower level is set to the maximum level of the fanouts
  template<typename TopoNtk>
  [[nodiscard]] std::pair<uint32_t, uint32_t> compute_level_span( const TopoNtk& topo_ntk, const std::vector<node<TopoNtk>>& nodes ) const noexcept
  {
    auto max_outgoing_level = topo_ntk.level( nodes.back() );
    topo_ntk.foreach_fanout( nodes.back(), [&topo_ntk, &max_outgoing_level]( const auto& fo )
                             {
                                const auto l = topo_ntk.level( fo );
                                if ( l > max_outgoing_level )
                                {
                                  max_outgoing_level = l;
                                } } );

    // the upper end of the span is the lowest level in the block
    const auto u = topo_ntk.level( nodes.front() );
    // the lower end of the span is the highest level in the block - 1 (by considering level gaps as dummy vertices of the block)
    const auto l = max_outgoing_level - 1;

    return { u, l };
  }

  // creates a single block from a list of nodes
  template<typename TopoNtk>
  [[nodiscard]] std::shared_ptr<block> create_block( const TopoNtk& topo_ntk, const node<TopoNtk>& n ) noexcept
  {
    std::vector<node<TopoNtk>> nodes{};
    gather_block_nodes( topo_ntk, n, nodes );

    auto b = std::make_shared<block>();
    b->nodes = nodes;
    b->span = compute_level_span( topo_ntk, nodes );

    std::for_each( b->nodes.cbegin(), b->nodes.cend(), [this, &b]( const auto& block_node )
                   { get_block[block_node] = b; } );

    return b;
  }

  // creates the initial block list
  void initialize_block_list() noexcept
  {
    ntk.incr_trav_id();

    topo_view topo_ntk{ ntk };
    topo_ntk.foreach_node( [this, &topo_ntk]( const auto& n )
                           {
                        // skip constants
                        if ( topo_ntk.is_constant( n ) )
                        {
                          return;
                        }
                        // skip nodes already visited
                        if ( topo_ntk.visited( n ) == topo_ntk.trav_id() )
                        {
                          return;
                        }

                        // create the block for node n
                        B.push_back( create_block( topo_ntk, n ) );
                        // store its current position within the block
                        B.back()->pi = B.size() - 1; } );
  }

  // don't deep-copy the block list to avoid unnecessary copies and invalidating the get_block associations
  [[nodiscard]] block_list copy_block_list_with_new_front( const std::shared_ptr<block>& new_front ) const noexcept
  {
    block_list B_prime{};
    B_prime.reserve( B.size() );

    B_prime.push_back( new_front );

    for ( const auto& b : B )
    {
      if ( b != new_front )
      {
        B_prime.push_back( b );
      }
    }

    return B_prime;
  }

  // moves the first block to index pos via multiple swaps
  void swap_front_to_position( block_list& b_list, const std::size_t pos ) const noexcept
  {
    for ( std::size_t i = 0u; i < pos; ++i )
    {
      assert( i < b_list.size() - 1 && "swap_front_to_position: pos out of bounds" );
      std::swap( b_list[i], b_list[i + 1] );
    }

    // TODO do we need to update the pi indices? probably not because sort_adjacencies will do that
  }

  void global_sifting() noexcept
  {
    // make a shallow copy of B for the subsequent iteration to ensure validity while B is being altered
    const auto B_copy = B;

    // for each sifting round
    for ( auto i = 0; i < params.rho; ++i )
    {
      for ( const auto& A : B_copy )
      {
        B = sifting_step( A );
      }
    }
  }

  block_list sifting_step( const std::shared_ptr<block>& A ) noexcept
  {
    // new block list has A in the front
    auto B_prime = copy_block_list_with_new_front( A );

    sort_adjacencies( B_prime );

    int32_t Chi = 0, Chi_star = 0; // current and best number of crossings
    std::size_t p_star = 0;        // best block position

    // start at 1 to skip A in B_prime
    for ( std::size_t p = 1; p < B_prime.size(); ++p )
    {
      Chi += sifting_swap( B_prime, A, B_prime[p] ); // these are consecutive blocks in B_prime

      // if a lower crossing number has been found, update best results
      if ( Chi < Chi_star )
      {
        Chi_star = Chi;
        p_star = p;
      }
    }

    // new block list ordered B'[1] ... B'[p*], A, B'[p* + 1] ... B'[ |B'| - 1 ]
    swap_front_to_position( B_prime, p_star );

    return B_prime;
  }

  // updates the N and I arrays
  void sort_adjacencies( const block_list& B_prime ) noexcept
  {
    // update pi and clear N and I
    for ( auto i = 0; i < B_prime.size(); ++i )
    {
      B_prime[i]->pi = i;

      // NOTE there is a hint in the paper that only updating the changed N and I improves runtime
      B_prime[i]->N_minus.clear();
      B_prime[i]->N_plus.clear();
      B_prime[i]->I_minus.clear();
      B_prime[i]->I_plus.clear();
    }

    // stores p values for each edge
    std::unordered_map<std::pair<node<Ntk>, node<Ntk>>, std::size_t, hash<std::pair<node<Ntk>, node<Ntk>>>> p{};

    // for each A in B'
    std::for_each( B_prime.cbegin(), B_prime.cend(), [this, &p]( const auto& A )
                   {
                     // update incoming adjacencies
                     const auto v = A->upper();
                     ntk.foreach_fanin( v, [this, &p, &A, &v]( const auto& fin )
                                        {
                                          const auto u = ntk.get_node( fin );
                                          const auto s = std::make_pair( u, v );
                                          const auto u_block = get_block[u];
                                          // add v to the next free position j of N^+(u)
                                          const auto j = u_block->N_plus.size();
                                          u_block->N_plus.push_back( v );

                                          if ( A->pi < u_block->pi ) // first traversal of segment s = (u, v)
                                          {
                                            p[s] = j;
                                          }
                                          else  // second traversal of segment s = (u, v)
                                          {
                                            u_block->I_plus[j] = p[s];
                                            get_block[v]->I_minus[p[s]] = j;
                                          } } );

                     // update outgoing adjacencies
                     const auto w = A->lower();
                     ntk.foreach_fanout( w, [this, &p, &A, &w]( const auto& x )
                                        {
                                          const auto s = std::make_pair( w, x );
                                          const auto x_block = get_block[x];
                                          // add w to the next free position j of N^-(x)
                                          const auto j = x_block->N_minus.size();
                                          x_block->N_minus.push_back( w );
                                          if ( A->pi < x_block->pi ) // first traversal of segment s = (w, x)
                                          {
                                            p[s] = j;
                                          }
                                          else  // second traversal of segment s = (w, x)
                                          {
                                            x_block->I_minus[j] = p[s];
                                            get_block[w]->I_plus[p[s]] = j;
                                          } } ); } );
  }

  enum direction
  {
    PLUS,
    MINUS
  };

  // returns the change in crossing count; A_block and B_block are consecutive blocks in b_list
  [[nodiscard]] int32_t sifting_swap( block_list& b_list, const std::shared_ptr<block>& A_block, const std::shared_ptr<block>& B_block ) noexcept
  {
    // set of level-direction pairs
    std::set<std::pair<uint32_t, direction>> L{};
    int32_t Delta{ 0 };

    // determine directions
    if ( const auto l = phi( A_block->upper() ); is_in_span( l, B_block ) )
    {
      L.insert( { l, direction::MINUS } );
    }
    if ( const auto l = phi( A_block->lower() ); is_in_span( l, B_block ) )
    {
      L.insert( { l, direction::PLUS } );
    }
    if ( const auto l = phi( B_block->upper() ); is_in_span( l, A_block ) )
    {
      L.insert( { l, direction::MINUS } );
    }
    if ( const auto l = phi( B_block->lower() ); is_in_span( l, A_block ) )
    {
      L.insert( { l, direction::PLUS } );
    }

    for ( const auto& [l, d] : L )
    {
      const auto a = get_node_by_level( A_block, l );
      const auto b = get_node_by_level( B_block, l );

      Delta += uswap( a, b, d );
      update_adjacencies( a, b, d );
    }

    // swap A_block and B_block in B_list
    std::swap( b_list[A_block->pi], b_list[B_block->pi] );

    ++( A_block->pi );
    --( B_block->pi );

    return Delta;
  }

  // determines the crossing count change for swapping a and b
  [[nodiscard]] int32_t uswap( const node<Ntk>& a, const node<Ntk>& b, const direction& d ) noexcept
  {
    const auto& x = d == direction::PLUS ? get_block[a]->N_plus : get_block[a]->N_minus;
    const auto& y = d == direction::PLUS ? get_block[b]->N_plus : get_block[b]->N_minus;

    int32_t c{ 0 };
    std::size_t i{ 0 }, j{ 0 };

    while ( i < x.size() && j < y.size() )
    {
      if ( get_block[x[i]]->pi < get_block[y[j]]->pi )
      {
        c = c + ( y.size() - j );
        i++;
      }
      else if ( get_block[x[i]]->pi > get_block[y[j]]->pi )
      {
        c = c - ( x.size() - i );
        j++;
      }
      else
      {
        c = c + ( y.size() - j ) - ( x.size() - i );
        i++;
        j++;
      }
    }

    return c;
  }

  void update_adjacencies( const node<Ntk>& a, const node<Ntk>& b, const direction& d ) noexcept
  {
    const auto& x = d == direction::PLUS ? get_block[a]->N_plus : get_block[a]->N_minus;
    const auto& y = d == direction::PLUS ? get_block[b]->N_plus : get_block[b]->N_minus;

    std::size_t i{ 0 }, j{ 0 };
    while ( i < x.size() && j < y.size() )
    {
      if ( get_block[x[i]]->pi < get_block[y[j]]->pi ) // x[i] is left of y[j]
      {
        i++;
      }
      else if ( get_block[x[i]]->pi > get_block[y[j]]->pi ) // x[i] is right of y[j]
      {
        j++;
      }
      else // x[i] is y[j]
      {
        const auto z = x[i]; // equal to y[j] because x[i] == y[j]

        // I^d(a)
        auto& a_I = ( d == direction::PLUS ? get_block[a]->I_plus : get_block[a]->I_minus );
        // I^d(b)
        auto& b_I = ( d == direction::PLUS ? get_block[b]->I_plus : get_block[b]->I_minus );

        // N^-d(z)
        auto& z_N = d == direction::PLUS ? get_block[z]->N_minus : get_block[z]->N_plus; // note the switched direction
        // I^-d(z)
        auto& z_I = d == direction::PLUS ? get_block[z]->I_minus : get_block[z]->I_plus; // note the switched direction

        // swap entries at positions I^d(a)[i] and I^d(b)[j] in N^-d(z) and in I^-d(z)
        std::swap( z_N[a_I[i]], z_N[b_I[j]] );
        std::swap( z_I[a_I[i]], z_I[b_I[j]] );

        a_I[i]++;
        b_I[j]--;

        i++;
        j++;
      }
    }
  }

  [[nodiscard]] uint32_t phi( const node<Ntk>& n ) const
  {
    return ntk.level( n );
  }

  // returns whether level is in the span of block b
  [[nodiscard]] bool is_in_span( const uint32_t level, const std::shared_ptr<block>& b ) const
  {
    return level >= b->span.first && level <= b->span.second;
  }

  // fetch a node from a block whose level is equal to the given level
  // if no such node exists, return the one with the smallest level difference
  [[nodiscard]] node<Ntk> get_node_by_level( const std::shared_ptr<block>& b, const uint32_t level ) const
  {
    node<Ntk> best_n{};
    uint32_t best_diff{ std::numeric_limits<uint32_t>::max() };

    for ( const auto& n : b->nodes )
    {
      if ( const auto diff = std::abs( static_cast<int32_t>( phi( n ) ) - static_cast<int32_t>( level ) ); diff < best_diff )
      {
        best_n = n;
        best_diff = diff;
      }
    }
    assert( best_diff != std::numeric_limits<uint32_t>::max() && "could not retrieve a node from an empty block" );

    return best_n;
  }
};

} // namespace detail

/**
 * @brief Optimizes the crossing count of a network by reordering its ranks.
 *
 * Based on "A Global k-Level Crossing Reduction Algorithm" by Christian Bachmaier, Franz J. Brandenburg, Wolfgang Brunner, and Ferdinand Hübner
 * in International Workshop on Algorithms and Computation (pp. 70-81).
 *
 * @tparam Ntk
 * @param rank_ntk
 * @param ps
 * @param pst
 */
template<typename Ntk>
void crossing_optimization( Ntk& rank_ntk, crossing_optimization_params const& ps = {}, crossing_optimization_stats* pst = nullptr ) noexcept
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node function" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin function" );
  static_assert( has_foreach_fanout_v<Ntk>, "Ntk does not implement the foreach_fanout function" );
  static_assert( has_fanin_size_v<Ntk>, "Ntk does not implement the fanin_size function" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size function" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node function" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant function" );
  static_assert( has_visited_v<Ntk>, "Ntk does not implement the visited function" );
  static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited function" );
  static_assert( has_trav_id_v<Ntk>, "Ntk does not implement the trav_id function" );
  static_assert( has_rank_position_v<Ntk>, "Ntk does not implement the rank_position function" );
  static_assert( has_level_v<Ntk>, "Ntk does not implement the level function" );

  crossing_optimization_stats st;
  detail::crossing_optimization_impl<Ntk> p{ rank_ntk, ps, st };
  p.run();

  if ( pst )
  {
    *pst = st;
  }
}

} // namespace mockturtle
