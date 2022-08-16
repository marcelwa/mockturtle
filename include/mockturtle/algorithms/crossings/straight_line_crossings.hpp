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
\file straight_line_crossings.hpp
\brief compute the straight-line crossing number of a network wrapped in a rank_view
\author Marcel Walter
*/

#pragma once

#include "../../traits.hpp"
#include "../../utils/stopwatch.hpp"

namespace mockturtle
{

struct straight_line_crossings_params
{
  int placeholder;
};

struct straight_line_crossings_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };
  /*! \brief Number of crossings */
  uint32_t num_crossings{ 0 };
};

namespace detail
{

template<typename Ntk>
class straight_line_crossings_impl
{
public:
  straight_line_crossings_impl( Ntk const& src, straight_line_crossings_params const& ps, straight_line_crossings_stats& st ) noexcept
      : ntk{ src },
        params{ ps },
        stats{ st }
  {
  }

  void run() noexcept
  {
    stopwatch t{ stats.time_total };

    foreach_edge( ntk, [this]( auto const& src1, const auto& tgt1 )
                  { foreach_edge( ntk, [this, &src1, &tgt1]( auto const& src2, const auto& tgt2 )
                                  {
                                    if ( is_straight_line_crossing( src1, tgt1, src2, tgt2 ) )
                                    {
                                      ++stats.num_crossings;
                                    } } ); } );

    // account for double counting of all crossings
    stats.num_crossings /= 2;
  }

private:
  /**
   * Network.
   */
  Ntk const& ntk;
  /**
   * Parameters.
   */
  straight_line_crossings_params const& params;
  /**
   * Statistics.
   */
  straight_line_crossings_stats& stats;

  template<typename Fun>
  void foreach_edge( Ntk const& ntk, Fun&& fun )
  {
    ntk.foreach_node( [&ntk, &fun]( auto const& tgt )
                      { ntk.foreach_fanin( tgt, [&ntk, &fun, &tgt]( auto const& fin )
                                           {
                                             if (auto const src = ntk.get_node( fin ); !ntk.is_constant(src))
                                             {
                                               fun( src, tgt );
                                             } } ); } );
  }

  struct node_pos
  {
    double x, y;
  };

  [[nodiscard]] node_pos get_node_pos( node<Ntk> const& n ) const noexcept
  {
    return { static_cast<double>( ntk.rank_position( n ) ), static_cast<double>( ntk.level( n ) ) };
  }
  /**
   * Computes whether two lines share a crossing point. The lines are defined by the endpoints of the given logic network nodes.
   *
   * This function uses Cramer's Rule to solve a linear equation system. The code is based on: https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
   *
   * @param src1 The source node of the first line.
   * @param tgt1 The target node of the first line.
   * @param src2 The source node of the second line.
   * @param tgt2 The target node of the second line.
   * @return True if the lines share a crossing point.
   */
  [[nodiscard]] bool is_straight_line_crossing( node<Ntk> const& src1, node<Ntk> const& tgt1, node<Ntk> const& src2, node<Ntk> const& tgt2 ) const noexcept
  {
    auto const p_src1 = get_node_pos( src1 );
    auto const p_tgt1 = get_node_pos( tgt1 );
    auto const p_src2 = get_node_pos( src2 );
    auto const p_tgt2 = get_node_pos( tgt2 );

    double const s1_x = p_tgt1.x - p_src1.x;
    double const s1_y = p_tgt1.y - p_src1.y;
    double const s2_x = p_tgt2.x - p_src2.x;
    double const s2_y = p_tgt2.y - p_src2.y;

    double const s = ( -s1_y * ( p_src1.x - p_src2.x ) + s1_x * ( p_src1.y - p_src2.y ) ) / ( -s2_x * s1_y + s1_x * s2_y );
    double const t = ( s2_x * ( p_src1.y - p_src2.y ) - s2_y * ( p_src1.x - p_src2.x ) ) / ( -s2_x * s1_y + s1_x * s2_y );

    // do not use >= / <= here to not consider lines that share the same endpoints as a crossing
    return ( s > 0 && s < 1 && t > 0 && t < 1 );
  }
};

} // namespace detail

template<typename Ntk>
void straight_line_crossings( Ntk const& rank_ntk, straight_line_crossings_params const& ps = {}, straight_line_crossings_stats* pst = nullptr ) noexcept
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node function" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin function" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node function" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant function" );
  static_assert( has_rank_position_v<Ntk>, "Ntk does not implement the rank_position function" );
  static_assert( has_level_v<Ntk>, "Ntk does not implement the level function" );

  straight_line_crossings_stats st;
  detail::straight_line_crossings_impl<Ntk> p{ rank_ntk, ps, st };
  p.run();

  if ( pst )
  {
    *pst = st;
  }
}

} // namespace mockturtle
