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
  \file spfd_resub.hpp
  \brief Simulation-guided resubstitution driven by SPFD statistical support selection.

  Same loop as `sim_resubstitution` (simulation + SAT validation with counter-example
  refinement); only the resynthesis engine differs.  Runtime knobs of the policy live in
  `spfd_global_params()`.
*/

#pragma once

#include "resyn_engines/spfd_resyn.hpp"
#include "sim_resub.hpp"

namespace mockturtle
{

template<class Ntk>
void spfd_sim_resubstitution( Ntk& ntk, resubstitution_params const& ps = {}, resubstitution_stats* pst = nullptr )
{
  static_assert( std::is_same_v<typename Ntk::base_type, aig_network> || std::is_same_v<typename Ntk::base_type, xag_network>,
                 "spfd_sim_resubstitution currently supports AIG and XAG" );

  using resub_view_t = fanout_view<depth_view<Ntk>>;
  depth_view<Ntk> depth_view{ ntk };
  resub_view_t resub_view{ depth_view };

  if constexpr ( std::is_same_v<typename Ntk::base_type, aig_network> )
  {
    using resyn_engine_t = spfd_resyn<kitty::partial_truth_table, aig_resyn_static_params_for_sim_resub<resub_view_t>>;
    using validator_t = circuit_validator<resub_view_t, bill::solvers::bsat2, false, true, false>;
    using resub_impl_t = typename detail::resubstitution_impl<resub_view_t, typename detail::simulation_based_resub_engine<resub_view_t, validator_t, resyn_engine_t>>;
    detail::sim_resubstitution_run<resub_view_t, resub_impl_t>( resub_view, ps, pst );
  }
  else
  {
    using resyn_engine_t = spfd_resyn<kitty::partial_truth_table, xag_resyn_static_params_for_sim_resub<resub_view_t>>;
    using validator_t = circuit_validator<resub_view_t, bill::solvers::bsat2, false, true, false>;
    using resub_impl_t = typename detail::resubstitution_impl<resub_view_t, typename detail::simulation_based_resub_engine<resub_view_t, validator_t, resyn_engine_t>>;
    detail::sim_resubstitution_run<resub_view_t, resub_impl_t>( resub_view, ps, pst );
  }
}

} /* namespace mockturtle */
