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
 \file crossing_graph.hpp
 \brief wrapper to use logic networks with the QuickCross library
 \author Marcel Walter
*/

#pragma once

#include "../../traits.hpp"

#include <cstdlib>

namespace mockturtle
{

template<typename Ntk>
class crossing_graph
{
public:
  explicit crossing_graph( Ntk const& ntk ) noexcept : network{ ntk }, num_vertices{ static_cast<int>( ntk.num_gates() + ntk.num_pis() + 1 ) }
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_num_gates_v<Ntk>, "Ntk does not implement the num_gates function" );
    static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node function" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin function" );
    static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index function" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node function" );
    static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant function" );

    init();
  }
  /**
   * The destructor frees the allocated edge_list memory.
   */
  ~crossing_graph()
  {
    free( edge_list );
  }
  /**
   * Returns the number of vertices in the graph.
   */
  [[nodiscard]] int get_num_vertices() const noexcept
  {
    return num_vertices;
  }
  /**
   * Returns the number of edges in the graph.
   */
  [[nodiscard]] int get_num_edges() const noexcept
  {
    return num_edges;
  }
  /**
   * Returns the edge list representation for the QuickCross library.
   *
   * @return the edge list representation for the QuickCross library.
   */
  [[nodiscard]] int* get_edge_list() const noexcept
  {
    return edge_list;
  }

private:
  /**
   * Store the associated network.
   */
  Ntk const& network;
  /**
   * Number of vertices in the crossing graph.
   */
  const int num_vertices;
  /**
   * Number of edges in the crossing graph.
   */
  int num_edges{ 0 };
  /**
   * An edge list representation of the graph.
   *
   * This "data structure" is used by the QuickCross library. It represents an unordered list of pairs of vertices.
   * To associate the vertices with each other into pairs, a simple index-based mapping is used. The first vertex of
   * each pair is located at position i with the second vertex at position i + M, where M is the number of all edges
   * (i.e., num_edges). Thus, the size of edge_list is '2 * M * sizeof(int)'.
   */
  int* edge_list;
  /**
   * Initialize the crossing graph.
   */
  void init() noexcept
  {
    // a lambda that applies a functor to each edge in a network
    const auto foreach_edge = []( auto const& ntk, auto&& fun ) -> void
    {
      ntk.foreach_node( [&ntk, &fun]( auto const& gate )
                        {
                          // no need to check for constants as they do not have fanins
                          ntk.foreach_fanin( gate, [&ntk, &fun, &gate]( auto const& fanin ) -> void
                                             {
                                               if (!ntk.is_constant(ntk.get_node(fanin)))
                                               {
                                                 fun( gate, fanin );
                                               }
                                             } ); } );
    };

    // determine the number of network edges
    foreach_edge( network, [this]( auto const& gate, auto const& fanin ) noexcept
                  { ++num_edges; } );

    // allocate memory for the edge list
    edge_list = static_cast<int*>( malloc( 2 * num_edges * sizeof( int ) ) );

    // fill the edge list
    foreach_edge( network, [this, i = 0u]( auto const& gate, auto const& fanin ) mutable noexcept
                  {
                    edge_list[i] = static_cast<int>(network.node_to_index(gate));
                    edge_list[i + num_edges] = static_cast<int>(network.node_to_index(network.get_node(fanin)));
                    ++i; } );
  }
};

} // namespace mockturtle