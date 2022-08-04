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

namespace quickcross
{
#include <quickcross/QC.c>
} // namespace quickcross

namespace mockturtle
{

struct crossing_number_params
{
  enum minimization_scheme : int
  {
    BIGFACE = 1,
    INCREMENTAL = 2,
    EXHAUSTIVE = 3
  };

  enum embedding_scheme : int
  {
    KAMADA_KAWAI_SPRING_MODEL = 1,
    CIRCLE_EMBEDDING = 2,
    PLANAR_EMBEDDING = 3,
    // NOTE: the following are not supported yet
    CROSSING_ORDER_LIST = 4,
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
   * If min_scheme == BIGFACE, this parameter determines the maximum search depth.
   */
  uint16_t bigface_depth{ 0 };
  /**
   * The initial graph embedding to start with.
   * KAMADA_KAWAI_SPRING_MODEL: The Kamada-Kawai spring embedding.
   * CIRCLE_EMBEDDING: Embedding onto a circle.
   * PLANAR_EMBEDDING: Iterative embedding upon a chordless cycle.
   */
  embedding_scheme initial_embedding{ PLANAR_EMBEDDING };
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

  // TODO coordinate list embedding
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
        qc_x_coordinates{ ps.initial_embedding == crossing_number_params::embedding_scheme::COORDINATE_LIST ? static_cast<long double*>( calloc( cg.get_num_vertices(), sizeof( long double ) ) ) : nullptr },
        qc_y_coordinates{ ps.initial_embedding == crossing_number_params::embedding_scheme::COORDINATE_LIST ? static_cast<long double*>( calloc( cg.get_num_vertices(), sizeof( long double ) ) ) : nullptr }
  {
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

    quickcross::Biconnected_Runner( cr_graph.get_edge_list(),                     // graph representation
                                    cr_graph.get_num_vertices(),                  // number of vertices N
                                    cr_graph.get_num_edges(),                     // number of edges M
                                    static_cast<int>( params.min_scheme ),        // minimization scheme
                                    static_cast<int>( params.initial_embedding ), // initial graph embedding
                                    static_cast<int>( params.bigface_depth ),     // max search depth for BIGFACE minimization
                                    static_cast<int>( params.random_seed ),       // random seed
                                    static_cast<int>( params.verbose ),           // be verbose
                                    0,                                            // output scheme, not being used in Biconnected_Runner
                                    0,                                            // stop iteration if a "sufficient" crossing number is found, not to be used
                                    nullptr,                                      // output file pointer, not being used in Biconnected_Runner
                                    qc_crossing_number,                           // stores the final crossing number
                                    &qc_crossing_label,
                                    qc_crossing_index,
                                    qc_x_coordinates,
                                    qc_y_coordinates );

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
