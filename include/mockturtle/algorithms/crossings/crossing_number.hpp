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
\file crossing_number.hpp
\brief use QuickCross to (approximately) compute the crossing number of a network wrapped in a crossing_graph
\author Marcel Walter
*/

#pragma once

#include "../../traits.hpp"
#include "../../utils/stopwatch.hpp"
#include "crossing_graph.hpp"

#include <cstdint>
#include <cstdlib>
#include <limits>

namespace quickcross
{
#include <quickcross/QC.h>
} // namespace quickcross

namespace mockturtle
{

struct crossing_number_params
{
  enum minimization_scheme : int
  {
    INCREMENTAL = 2,
    EXHAUSTIVE = 3
  };

  enum embedding_scheme : int
  {
    PLANAR_EMBEDDING = 3,
    COORDINATE_LIST = 5
  };
  /**
   * The crossing minimization technique to use.
   * BIGFACE: Start with the graph face of the most sides and iteratively deepen from there, rearranging vertices to find crossing number improvements. Breaks as soon as an improvement is found.
   * INCREMENTAL: Standard QuickCross, breaks until an improvement has been found.
   * EXHAUSTIVE: Standard QuickCross, exhaustively looks at all improvements and settles for the best.
   */
  minimization_scheme min_scheme{ INCREMENTAL };
  /**
   * The initial graph embedding to start with.
   * PLANAR_EMBEDDING: Iterative embedding upon a chordless cycle.
   * COORDINATE_LIST: Straight-line embedding given a coordinate list (extracted from rank_view).
   */
  embedding_scheme initial_embedding{ PLANAR_EMBEDDING };
  /**
   * Execute a single run or perform optimization.
   */
  bool single_run{ true };
  /**
   * Specifies a seed for embedding randomization. If random_seed == 0, no randomization is performed.
   */
  uint32_t random_seed{ 0 };
  /**
   * Be verbose.
   */
  bool verbose{ false };
};

struct crossing_number_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };
  /*! \brief Number of crossings */
  uint32_t num_crossings{ 0 };
};

namespace detail
{

template<typename Ntk>
class crossing_number_impl
{
public:
  crossing_number_impl( crossing_graph<Ntk> const& cg, crossing_number_params const& ps, crossing_number_stats& st ) noexcept
      : cr_graph{ cg },
        params{ ps },
        stats{ st },
        qc_crossing_label{ static_cast<int*>( calloc( cg.get_num_edges(), sizeof( int ) ) ) },
        qc_crossing_index{ static_cast<int*>( calloc( cg.get_num_edges() + 1, sizeof( int ) ) ) },
        qc_x_coordinates{ ( ps.initial_embedding == crossing_number_params::embedding_scheme::COORDINATE_LIST && has_level_v<Ntk> && has_rank_position_v<Ntk> ) ? static_cast<long double*>( calloc( cg.get_num_vertices(), sizeof( long double ) ) ) : nullptr },
        qc_y_coordinates{ ( ps.initial_embedding == crossing_number_params::embedding_scheme::COORDINATE_LIST && has_level_v<Ntk> && has_rank_position_v<Ntk> ) ? static_cast<long double*>( calloc( cg.get_num_vertices(), sizeof( long double ) ) ) : nullptr }
  {
    if constexpr ( has_level_v<Ntk> && has_rank_position_v<Ntk> )
    {
      if ( ps.initial_embedding == crossing_number_params::embedding_scheme::COORDINATE_LIST )
      {
        cr_graph.ntk().foreach_node( [this]( auto const& n )
                                     {
                                        auto const idx = cr_graph.ntk().node_to_index( n );

                                        if ( n == cr_graph.ntk().get_node( cr_graph.ntk().get_constant( false ) ) )
                                        {
                                          qc_x_coordinates[idx] = static_cast<long double>( cr_graph.ntk().depth() + 2 );
                                          qc_y_coordinates[idx] = static_cast<long double>( cr_graph.ntk().width() + 2 );
                                        }
                                        else if ( n == cr_graph.ntk().get_node( cr_graph.ntk().get_constant( true ) ) )
                                        {
                                          qc_x_coordinates[idx] = static_cast<long double>( cr_graph.ntk().depth() + 3 );
                                          qc_y_coordinates[idx] = static_cast<long double>( cr_graph.ntk().width() + 3 );
                                        }
                                        else
                                        {
                                          qc_x_coordinates[idx] = static_cast<long double>( cr_graph.ntk().rank_position( n ) + 1 );
                                          qc_y_coordinates[idx] = static_cast<long double>( cr_graph.ntk().level( n ) + 1);
                                        } } );
      }
    }
  }

  ~crossing_number_impl()
  {
    free( qc_crossing_label );
    free( qc_crossing_index );
    free( qc_x_coordinates );
    free( qc_y_coordinates );
  }

  void run() noexcept
  {
    stopwatch t{ stats.time_total };

    quickcross::Biconnected_Runner( cr_graph.get_edge_list(),                                // graph representation
                                    cr_graph.get_num_vertices(),                             // number of vertices N
                                    cr_graph.get_num_edges(),                                // number of edges M
                                    static_cast<int>( params.min_scheme ),                   // minimization scheme
                                    static_cast<int>( params.initial_embedding ),            // initial graph embedding
                                    0,                                                       // max search depth for BIGFACE minimization, not being used by us
                                    static_cast<int>( params.random_seed ),                  // random seed
                                    static_cast<int>( params.verbose ),                      // be verbose
                                    0,                                                       // output scheme, not being used in Biconnected_Runner
                                    params.single_run ? std::numeric_limits<int>::max() : 0, // stop iteration if a "sufficient" crossing number is found
                                    nullptr,                                                 // output file pointer, not being used in Biconnected_Runner
                                    qc_crossing_number,                                      // stores the final crossing number
                                    &qc_crossing_label,                                      // stores the final crossing labels
                                    qc_crossing_index,                                       // stores the final crossing indices
                                    qc_x_coordinates,                                        // list of x-coordinates of all vertices, only used if initial_embedding == COORDINATE_LIST
                                    qc_y_coordinates                                         // list of y-coordinates of all vertices, only used if initial_embedding == COORDINATE_LIST
    );

    // extract the crossing number from the algorithm output
    stats.num_crossings = static_cast<uint32_t>( qc_crossing_number[0] );
  }

private:
  /**
   * The network wrapper whose crossing number is computed.
   */
  crossing_graph<Ntk> const& cr_graph;
  /**
   * Parameters.
   */
  crossing_number_params const& params;
  /**
   * Statistics.
   */
  crossing_number_stats& stats;
  /**
   * Host the crossing number after the QuickCross run.
   */
  int qc_crossing_number[1] = { -1 };
  /**
   * Hosts the crossing labels after the QuickCross run.
   */
  int* qc_crossing_label;
  /**
   * Hosts the crossing indices after the QuickCross run.
   */
  int* qc_crossing_index;
  /**
   * List of x-coordinates for each vertex. Only used for the COORDINATE_LIST embedding.
   */
  long double* qc_x_coordinates;
  /**
   * List of y-coordinates for each vertex. Only used for the COORDINATE_LIST embedding.
   */
  long double* qc_y_coordinates;
};

} // namespace detail

template<typename Ntk>
void crossing_number( crossing_graph<Ntk> const& cg, crossing_number_params const& ps = {}, crossing_number_stats* pst = nullptr ) noexcept
{
  crossing_number_stats st;
  detail::crossing_number_impl<Ntk> p{ cg, ps, st };
  p.run();

  if ( pst )
  {
    *pst = st;
  }
}

} // namespace mockturtle
