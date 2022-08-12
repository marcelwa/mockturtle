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
};

} // namespace detail

template<typename Ntk>
void straight_line_crossings( Ntk const& rank_ntk, straight_line_crossings_params const& ps = {}, straight_line_crossings_stats* pst = nullptr ) noexcept
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_rank_position_v<Ntk>, "Ntk does not implement the rank_position function" );
  static_assert( has_at_rank_position_v<Ntk>, "Ntk does not implement the at_rank_position function" );

  straight_line_crossings_stats st;
  detail::straight_line_crossings_impl<Ntk> p{ rank_ntk, ps, st };
  p.run();

  if ( pst )
  {
    *pst = st;
  }
}

} // namespace mockturtle
