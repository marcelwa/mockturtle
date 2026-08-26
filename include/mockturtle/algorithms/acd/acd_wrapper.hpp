/* Thin C++ interface to the vendored Ashenhurst-Curtis decomposer.
 *
 * Mirrors ABC's src/map/if/acd/ac_wrapper.cpp so that the two implementations
 * ask the decomposer exactly the same questions.  Header-only.
 */
#pragma once

#include <cstdint>

#include "acd_namespace.hpp"
#include "ac_decomposition.hpp"

namespace mockturtle
{

namespace acd_iface
{

using word = abcacd::word;

/*! \brief Feasibility + cost of an AC decomposition of `ptt` into `lut_size`-LUTs.
 *
 * \param ptt      truth table of the cut, `num_vars` variables, ABC word layout
 * \param num_vars support size of the cut (must exceed `lut_size` to be useful)
 * \param lut_size target LUT size of the decomposed structure
 * \param cost     out: number of LUTs of the structure
 * \return number of levels of the structure, or -1 if infeasible
 *
 * The delay profile is always empty: no leaf is declared late-arriving, so the
 * decomposer is free to place any variable in the free set and the answer does
 * not depend on the mapper's current arrival times.  This makes the cost of a
 * cut a function of the cut alone, which is what an area-oriented comparison
 * wants; see the method chapter for why we deviate from ABC here.
 */
inline int evaluate( word* ptt, uint32_t num_vars, uint32_t lut_size, uint32_t* cost )
{
  using namespace abcacd::acd;

  ac_decomposition_params ps;
  ps.lut_size = lut_size;
  ps.use_first = false;
  ps.try_no_late_arrival = false;
  ac_decomposition_stats st;

  ac_decomposition_impl acd( num_vars, ps, &st );
  unsigned delay_profile = 0u;
  int const levels = acd.run( ptt, delay_profile );

  if ( levels < 0 )
  {
    *cost = 0;
    return -1;
  }

  *cost = st.num_luts;
  return levels;
}

/*! \brief Computes the decomposition and writes it into ABC's `decompArray` format.
 *
 * Layout (identical to ABC's, see Abc_DecRecordToHop):
 *   [0]      total number of bytes written
 *   [1]      number of LUTs, last one is the root
 *   then per LUT: number of fanins, then that many fanin indices
 *   (index < num_vars: cut leaf; else index - num_vars: earlier LUT of this
 *   structure), then the truth table, little-endian, one byte at a time.
 *
 * \return 0 on success, -1 if no decomposition exists.
 */
inline int decompose( word* ptt, uint32_t num_vars, uint32_t lut_size, unsigned char* decomp_array )
{
  using namespace abcacd::acd;

  ac_decomposition_params ps;
  ps.lut_size = lut_size;
  ps.use_first = true;
  ac_decomposition_stats st;

  ac_decomposition_impl acd( num_vars, ps, &st );
  unsigned delay_profile = 0u;
  if ( acd.run( ptt, delay_profile ) < 0 )
    return -1;
  if ( acd.compute_decomposition() < 0 )
    return -1;

  acd.get_decomposition( decomp_array );
  return 0;
}

} /* namespace acd_iface */

} /* namespace mockturtle */
