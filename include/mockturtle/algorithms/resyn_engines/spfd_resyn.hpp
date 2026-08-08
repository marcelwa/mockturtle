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
  \file spfd_resyn.hpp
  \brief Resynthesis with SPFD / information-graph statistical support selection.

  Implements the support-selection policy of

    A. Costamagna, A. Tempia Calvino, A. Mishchenko, G. De Micheli,
    "Enhanced Resubstitution for Logic Optimization" / "Area-Oriented Resynthesis
    with Information Graphs" (ISCAS'24, IWLS'24),

  on top of mockturtle's existing simulation-guided resubstitution loop.

  The engine is a *decorator*: it selects a small divisor support `C` by covering the
  information graph (SPFD) of the target, then hands `C` -- and only `C` -- to an inner
  `xag_resyn_decompose` for the actual dependency-circuit construction.  Because supports
  are *sampled* rather than taken greedily, `num_supports` independent supports give the
  inner engine `num_supports` structurally different chances at the same pivot.  Only the
  smallest resulting dependency circuit is returned, so this engine is never worse than
  the plain engine on a single call, only slower.

  \author agentic-synthesis (Costamagna et al.'s policy)
*/

#pragma once

#include "../../utils/index_list/index_list.hpp"
#include "../../utils/stopwatch.hpp"
#include "xag_resyn.hpp"

#include <fmt/format.h>
#include <kitty/kitty.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <vector>

namespace mockturtle
{

/*! \brief Runtime parameters of `spfd_resyn`.
 *
 * The resynthesis-engine interface of `resubstitution.hpp` constructs the engine with a
 * stats reference only, so runtime knobs are read from this process-global block.  Set it
 * once before calling `sim_resubstitution`.
 */
struct spfd_resyn_params
{
  /*! \brief `K` -- maximum support size explored by the covering process. */
  uint32_t max_support{ 7u };

  /*! \brief `S` -- number of supports sampled per resynthesis call. 0 disables SPFD. */
  uint32_t num_supports{ 10u };

  /*! \brief `beta` -- inverse temperature of the Boltzmann policy on the *normalised*
   *  remaining-edge count. A negative value means pure greedy selection (beta -> inf). */
  double beta{ 5.0 };

  /*! \brief Cap on the number of divisors entering the covering process. */
  uint32_t max_divisors{ 150u };

  /*! \brief Random seed. */
  uint32_t seed{ 1u };

  /*! \brief Only sample supports when the plain engine found nothing (cheaper). */
  bool only_on_fail{ false };

  /*! \brief Diagnostic: when a valid support does not fit the budget, retry it with an
   *  unlimited budget. Distinguishes "the budget is the binding constraint" from "the
   *  covering process is returning supports that are not actually valid". */
  bool diagnose{ false };
};

inline spfd_resyn_params& spfd_global_params()
{
  static spfd_resyn_params ps;
  return ps;
}

struct spfd_resyn_stats
{
  /*! \brief Time spent in the covering / support-sampling process. */
  stopwatch<>::duration time_cover{ 0 };

  /*! \brief Number of resynthesis calls. */
  uint32_t num_calls{ 0 };

  /*! \brief Calls in which the plain (full-divisor) engine succeeded. */
  uint32_t num_base_success{ 0 };

  /*! \brief Calls whose insertion budget was 0, i.e. the pivot's MFFC has one node so no
   *  support of size > 1 could ever be paid for. Diagnostic for why SPFD does or does not
   *  get a chance in an AIG resubstitution framework. */
  uint32_t num_calls_no_budget{ 0 };

  /*! \brief Sum of the effective support caps `k_eff` over all sampled supports. */
  uint64_t sum_k_eff{ 0 };

  /*! \brief Valid supports for which the inner engine could not fit the dependency
   *  circuit into the budget. */
  uint32_t num_valid_but_unfit{ 0 };

  /*! \brief Of those, how many are realisable at all (diagnostic mode only). */
  uint32_t num_unfit_realisable{ 0 };
  uint64_t sum_unfit_size{ 0 };
  uint64_t sum_unfit_budget{ 0 };

  /*! \brief Supports sampled. */
  uint32_t num_supports_sampled{ 0 };

  /*! \brief Supports that satisfied the dependency theorem on the signatures. */
  uint32_t num_supports_valid{ 0 };

  /*! \brief Calls where SPFD found a solution and the plain engine found none. */
  uint32_t num_spfd_only{ 0 };

  /*! \brief Calls where SPFD found a *strictly smaller* solution than the plain engine. */
  uint32_t num_spfd_smaller{ 0 };

  /*! \brief Total gates saved by SPFD relative to the plain engine (on smaller-solution calls). */
  uint32_t num_gates_saved{ 0 };

  xag_resyn_stats inner_st;

  void report() const
  {
    fmt::print( "[i]         <spfd_resyn>\n" );
    fmt::print( "[i]             #calls          : {:>8d}\n", num_calls );
    fmt::print( "[i]             #base success   : {:>8d}\n", num_base_success );
    fmt::print( "[i]             #zero budget    : {:>8d} (MFFC = 1, no support > 1 payable)\n", num_calls_no_budget );
    fmt::print( "[i]             #supports       : {:>8d} ({} valid, {} unfit, mean k_eff {:.2f})\n",
                num_supports_sampled, num_supports_valid, num_valid_but_unfit,
                num_supports_sampled ? double( sum_k_eff ) / double( num_supports_sampled ) : 0.0 );
    if ( num_unfit_realisable )
    {
      fmt::print( "[i]             unfit-but-realisable: {} (mean {:.2f} gates vs budget {:.2f})\n",
                  num_unfit_realisable, double( sum_unfit_size ) / double( num_unfit_realisable ),
                  double( sum_unfit_budget ) / double( num_unfit_realisable ) );
    }
    fmt::print( "[i]             #spfd-only wins : {:>8d}\n", num_spfd_only );
    fmt::print( "[i]             #spfd-smaller   : {:>8d} (-{} gates)\n", num_spfd_smaller, num_gates_saved );
    fmt::print( "[i]             covering        : {:>8.2f} secs\n", to_seconds( time_cover ) );
    inner_st.report();
  }
};

namespace detail
{

/*! \brief |a & b| without materialising `a & b`. */
template<class TT>
inline uint32_t spfd_count_intersection( TT const& a, TT const& b )
{
  uint32_t c = 0u;
  auto ia = a.begin();
  auto ib = b.begin();
  for ( ; ia != a.end(); ++ia, ++ib )
  {
    c += static_cast<uint32_t>( __builtin_popcountll( *ia & *ib ) );
  }
  return c;
}

} /* namespace detail */

/*! \brief Resynthesis engine with SPFD statistical support selection.
 *
 * ## The information graph
 *
 * Let `x` be the target and `x~` its `p`-bit simulation signature.  The SPFD (information
 * graph) of `x` under care set `c` is the complete bipartite graph joining every care
 * minterm of the on-set to every care minterm of the off-set; it therefore has
 * `|ON| * |OFF|` edges.  The dependency theorem states that a function `g` with
 * `x = g(C)` exists iff every such pair is distinguished by some divisor in `C`.
 *
 * ## The covering state
 *
 * Covering by a divisor `d` splits every current partition into the sub-partitions where
 * `d = 0` and where `d = 1`; edges *inside* a sub-partition remain uncovered.  Storing the
 * partition as the pair `(on_i, off_i)` of signature masks, the number of edges still
 * uncovered after additionally covering with `d` is
 *
 *   H(state, d) = sum_i  |on_i & d| * |off_i & d|
 *               +        (|on_i| - |on_i & d|) * (|off_i| - |off_i & d|),
 *
 * i.e. two popcounts per partition per divisor.  Partitions in which the target is
 * constant contribute zero edges forever and are dropped, which is what keeps the state
 * from growing as `2^t`.
 *
 * ## The policy
 *
 * With `H_cur = H(state)` and `H_min = min_d H(state, d)`, the normalised cost is
 *
 *   Hn(d) = (H(state,d) - H_min) / (H_cur - H_min)  in [0,1],
 *
 * and the next divisor is drawn with `p(d) ~ exp(-beta * Hn(d))`.  `beta -> inf` recovers
 * greedy support selection; `beta = 0` is uniform.
 */
template<class TT, class static_params = xag_resyn_static_params_default<TT>>
class spfd_resyn
{
public:
  using stats = spfd_resyn_stats;
  using index_list_t = large_xag_index_list;
  using truth_table_t = TT;
  using inner_engine_t = xag_resyn_decompose<TT, static_params>;

private:
  struct cover_part
  {
    TT on, off;
    uint32_t non, noff;
  };

public:
  explicit spfd_resyn( stats& st ) noexcept
      : st( st ), inner( st.inner_st ), ps( spfd_global_params() ), rng( spfd_global_params().seed )
  {
    static_assert( !static_params::copy_tts, "spfd_resyn requires copy_tts = false" );
  }

  /*! \brief Perform resynthesis.
   *
   * Same interface as `xag_resyn_decompose::operator()`.  Input literals of the returned
   * index list refer to the *full* divisor range `[begin, end)`, as the resubstitution
   * framework requires.
   */
  template<class iterator_type,
           bool enabled = static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end,
                                          typename static_params::truth_table_storage_type const& tts,
                                          uint32_t max_size = std::numeric_limits<uint32_t>::max() )
  {
    ++st.num_calls;
    ptts = &tts;
    divs.clear();
    for ( auto it = begin; it != end; ++it )
    {
      divs.push_back( *it );
    }

    /* 1. the plain engine on the full divisor set -- this is exactly the control arm */
    std::optional<index_list_t> best = inner( target, care, begin, end, tts, max_size );
    uint32_t best_size = best ? static_cast<uint32_t>( best->num_gates() ) : max_size + 1u;
    if ( best )
    {
      ++st.num_base_success;
    }

    if ( max_size == 0u )
    {
      /* Only a 0-resub is payable; the plain engine already enumerates those exhaustively
         and any support of size > 1 is unrealisable. Nothing for SPFD to do. */
      ++st.num_calls_no_budget;
      return best;
    }
    if ( ps.num_supports == 0u || best_size == 0u || divs.empty() )
    {
      return best;
    }
    if ( ps.only_on_fail && best )
    {
      return best;
    }

    /* 2. SPFD-sampled supports */
    on_set = target & care;
    off_set = ~target & care;
    if ( kitty::count_ones( on_set ) == 0u || kitty::count_ones( off_set ) == 0u )
    {
      return best; /* constant under the care set; the plain engine already handles it */
    }

    /* candidate divisors: the last `max_divisors` entries, which the collector orders
       closest to the pivot */
    cand.clear();
    uint32_t const n_divs = static_cast<uint32_t>( divs.size() );
    uint32_t const first = n_divs > ps.max_divisors ? n_divs - ps.max_divisors : 0u;
    for ( uint32_t i = first; i < n_divs; ++i )
    {
      cand.push_back( i );
    }

    bool const had_base = best.has_value();
    uint32_t const base_size = best_size;

    for ( uint32_t s = 0u; s < ps.num_supports; ++s )
    {
      if ( best_size == 0u )
      {
        break;
      }
      /* A dependency circuit of `g` two-input gates has at most `g + 1` distinct inputs, so
         a support larger than `budget + 1` can never be realised within the budget.  This
         caps the covering process at the useful support size instead of always at K. */
      uint32_t const budget = best ? best_size - 1u : max_size;
      uint32_t const k_eff = std::min<uint32_t>( ps.max_support, budget == std::numeric_limits<uint32_t>::max() ? ps.max_support : budget + 1u );
      if ( k_eff == 0u )
      {
        break;
      }

      ++st.num_supports_sampled;
      st.sum_k_eff += k_eff;
      bool const ok = call_with_stopwatch( st.time_cover, [&]() { return sample_support( k_eff ); } );
      if ( !ok )
      {
        continue;
      }
      ++st.num_supports_valid;

      supp_nodes.clear();
      for ( auto const i : support )
      {
        supp_nodes.push_back( divs[i] );
      }

      auto const res = inner( target, care, supp_nodes.begin(), supp_nodes.end(), tts, budget );
      if ( res )
      {
        best = remap( *res, n_divs );
        best_size = static_cast<uint32_t>( best->num_gates() );
      }
      else
      {
        ++st.num_valid_but_unfit;
        if ( ps.diagnose )
        {
          auto const free_res = inner( target, care, supp_nodes.begin(), supp_nodes.end(), tts, 64u );
          if ( free_res )
          {
            ++st.num_unfit_realisable;
            st.sum_unfit_size += static_cast<uint32_t>( free_res->num_gates() );
            st.sum_unfit_budget += budget;
          }
        }
      }
    }

    if ( best && !had_base )
    {
      ++st.num_spfd_only;
    }
    else if ( best && had_base && best_size < base_size )
    {
      ++st.num_spfd_smaller;
      st.num_gates_saved += base_size - best_size;
    }

    return best;
  }

private:
  /*! \brief Sample one support by the Boltzmann covering policy.
   *
   * Returns true and fills `support` (indices into `divs`) iff the sampled support covers
   * the whole information graph, i.e. satisfies the dependency theorem on the signatures.
   */
  bool sample_support( uint32_t k_eff )
  {
    support.clear();
    parts.clear();
    parts.push_back( cover_part{ on_set, off_set, kitty::count_ones( on_set ), kitty::count_ones( off_set ) } );
    uint64_t h_cur = static_cast<uint64_t>( parts[0].non ) * static_cast<uint64_t>( parts[0].noff );

    while ( support.size() < k_eff )
    {
      /* remaining-edge count for every candidate divisor */
      costs.clear();
      costs.reserve( cand.size() );
      uint64_t h_min = std::numeric_limits<uint64_t>::max();
      for ( auto const i : cand )
      {
        if ( std::find( support.begin(), support.end(), i ) != support.end() )
        {
          costs.push_back( std::numeric_limits<uint64_t>::max() );
          continue;
        }
        uint64_t const h = h_after( ( *ptts )[divs[i]] );
        costs.push_back( h );
        h_min = std::min( h_min, h );
      }

      if ( h_min >= h_cur )
      {
        return false; /* no divisor makes progress */
      }

      /* Boltzmann weights on the normalised cost */
      uint32_t pick = 0u;
      if ( ps.beta < 0.0 )
      {
        /* greedy, uniform tie-break */
        ties.clear();
        for ( uint32_t j = 0u; j < costs.size(); ++j )
        {
          if ( costs[j] == h_min )
          {
            ties.push_back( j );
          }
        }
        pick = ties[std::uniform_int_distribution<uint32_t>( 0u, static_cast<uint32_t>( ties.size() ) - 1u )( rng )];
      }
      else
      {
        double const denom = static_cast<double>( h_cur - h_min );
        double total = 0.0;
        weights.clear();
        weights.reserve( costs.size() );
        for ( auto const c : costs )
        {
          double w = 0.0;
          if ( c != std::numeric_limits<uint64_t>::max() && c < h_cur )
          {
            double const hn = static_cast<double>( c - h_min ) / denom;
            w = std::exp( -ps.beta * hn );
          }
          weights.push_back( w );
          total += w;
        }
        if ( total <= 0.0 )
        {
          return false;
        }
        double r = std::uniform_real_distribution<double>( 0.0, total )( rng );
        pick = static_cast<uint32_t>( weights.size() ) - 1u;
        for ( uint32_t j = 0u; j < weights.size(); ++j )
        {
          r -= weights[j];
          if ( r <= 0.0 )
          {
            pick = j;
            break;
          }
        }
      }

      uint32_t const chosen = cand[pick];
      support.push_back( chosen );
      h_cur = apply_cover( ( *ptts )[divs[chosen]] );
      if ( h_cur == 0u )
      {
        return true;
      }
    }
    return false;
  }

  /*! \brief Remaining edges after additionally covering the current state with `d`. */
  uint64_t h_after( TT const& d ) const
  {
    uint64_t h = 0u;
    for ( auto const& p : parts )
    {
      uint32_t const c1 = detail::spfd_count_intersection( p.on, d );
      uint32_t const c0 = p.non - c1;
      uint32_t const k1 = detail::spfd_count_intersection( p.off, d );
      uint32_t const k0 = p.noff - k1;
      h += static_cast<uint64_t>( c1 ) * k1 + static_cast<uint64_t>( c0 ) * k0;
    }
    return h;
  }

  /*! \brief Commit a divisor into the covering state; returns the new edge count. */
  uint64_t apply_cover( TT const& d )
  {
    next_parts.clear();
    uint64_t h = 0u;
    for ( auto const& p : parts )
    {
      TT on1 = p.on & d;
      TT off1 = p.off & d;
      uint32_t const c1 = kitty::count_ones( on1 );
      uint32_t const k1 = kitty::count_ones( off1 );
      if ( c1 != 0u && k1 != 0u )
      {
        h += static_cast<uint64_t>( c1 ) * k1;
        next_parts.push_back( cover_part{ std::move( on1 ), std::move( off1 ), c1, k1 } );
      }
      TT on0 = p.on & ~d;
      TT off0 = p.off & ~d;
      uint32_t const c0 = p.non - c1;
      uint32_t const k0 = p.noff - k1;
      if ( c0 != 0u && k0 != 0u )
      {
        h += static_cast<uint64_t>( c0 ) * k0;
        next_parts.push_back( cover_part{ std::move( on0 ), std::move( off0 ), c0, k0 } );
      }
    }
    parts.swap( next_parts );
    return h;
  }

  /*! \brief Re-index a dependency circuit built over `support` onto the full divisor set. */
  index_list_t remap( index_list_t const& in, uint32_t num_divs ) const
  {
    index_list_t out;
    out.add_inputs( num_divs );
    std::vector<uint32_t> gate_lit;
    uint32_t const n_in = static_cast<uint32_t>( in.num_pis() );

    auto map_lit = [&]( uint32_t lit ) -> uint32_t {
      uint32_t const idx = lit >> 1;
      uint32_t const comp = lit & 1u;
      if ( idx == 0u )
      {
        return comp; /* constant */
      }
      if ( idx <= n_in )
      {
        return ( ( support[idx - 1u] + 1u ) << 1 ) | comp;
      }
      return gate_lit[idx - n_in - 1u] | comp;
    };

    in.foreach_gate( [&]( uint32_t l0, uint32_t l1 ) {
      uint32_t const a = map_lit( l0 );
      uint32_t const b = map_lit( l1 );
      gate_lit.push_back( in.is_and( l0, l1 ) ? out.add_and( a, b ) : out.add_xor( a, b ) );
    } );

    uint32_t po = 0u;
    in.foreach_po( [&]( uint32_t l ) { po = map_lit( l ); } );
    out.add_output( po );
    return out;
  }

private:
  stats& st;
  inner_engine_t inner;
  spfd_resyn_params const& ps;
  std::mt19937 rng;

  const typename static_params::truth_table_storage_type* ptts{ nullptr };
  std::vector<typename static_params::node_type> divs;
  std::vector<typename static_params::node_type> supp_nodes;
  std::vector<uint32_t> cand;
  std::vector<uint32_t> support;
  std::vector<uint64_t> costs;
  std::vector<double> weights;
  std::vector<uint32_t> ties;
  std::vector<cover_part> parts, next_parts;
  TT on_set, off_set;
}; /* spfd_resyn */

} /* namespace mockturtle */
