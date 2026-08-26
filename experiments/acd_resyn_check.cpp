/* acd_resyn_check: the correctness gate for `node_resynthesis/acd.hpp`.
 *
 * No QoR number from the ACD resynthesis operator may be quoted before this passes.  It
 * follows the L3v protocol, which caught a vacuous "0 failures" in this project once
 * before (a `verified` flag initialised to 1): every emitted replacement is checked by
 * *exhaustive* 2^N simulation against the function it was asked to implement, and a
 * deliberate sabotage control must fire.  If the sabotage does not fire the checker is
 * broken and the run exits non-zero regardless of how many honest cases passed.
 *
 * Two families of functions are exercised.  Uniformly random truth tables are the worst
 * case for a decomposer and are what stresses the decompArray walk; but they are also
 * almost never AC-decomposable above 8 variables, so on their own they would leave the
 * 9-11 variable cascades -- exactly the ones a real circuit produces -- untested.  So the
 * second family is real: the cut functions of an actual AIG, which is what the operator
 * will be handed in anger.
 *
 * Usage: acd_resyn_check [--seed=N] [--reps=N] [--aig=PATH] [--verbose]
 */

#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include <lorina/aiger.hpp>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/algorithms/node_resynthesis/acd.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>

namespace mockturtle
{
/* mockturtle's cut priority list is smallest-first, which starves this arm of the wide
 * cuts it exists to test (see 2026-08-26-acd-resyn-census).  Reverse it. */
struct check_wide_cut
{
};
template<bool ComputeTruth>
bool operator<( cut_type<ComputeTruth, check_wide_cut> const& c1,
                cut_type<ComputeTruth, check_wide_cut> const& c2 )
{
  return c1.size() > c2.size();
}
} /* namespace mockturtle */

using namespace mockturtle;

namespace
{

struct outcome
{
  uint32_t calls{ 0 };
  uint32_t emitted{ 0 };   /* the functor invoked its callback                    */
  uint32_t accepted{ 0 };  /* ACD produced the cascade (as opposed to a fallback) */
  uint32_t mismatch{ 0 };  /* emitted, but not equal to the requested function    */
};

/*! \brief Runs one resynthesis call in a fresh network and checks it exhaustively.
 *
 * \param sabotage when true the *checker's* reference function is left alone but the
 *        emitted signal is complemented, so a correct checker must report a mismatch.
 */
template<class Ntk, class ResynFn>
bool check_one( ResynFn& resyn, kitty::dynamic_truth_table const& target,
                bool sabotage, bool* emitted, bool* by_acd,
                acd_resynthesis_stats const& before, acd_resynthesis_stats const& probe )
{
  uint32_t const n = target.num_vars();

  Ntk ntk;
  std::vector<signal<Ntk>> pis;
  for ( uint32_t i = 0; i < n; ++i )
    pis.push_back( ntk.create_pi() );

  bool got = false;
  signal<Ntk> out = ntk.get_constant( false );
  resyn( ntk, target, pis.begin(), pis.end(), [&]( signal<Ntk> const& s ) {
    if ( !got )
    {
      out = s;
      got = true;
    }
    return true;
  } );

  *emitted = got;
  *by_acd = ( probe.num_accepted != before.num_accepted );
  if ( !got )
    return true; /* declining is allowed; emitting something wrong is not */

  if ( sabotage )
    out = ntk.create_not( out );

  ntk.create_po( out );

  default_simulator<kitty::dynamic_truth_table> sim( static_cast<unsigned>( n ) );
  auto const tts = simulate<kitty::dynamic_truth_table>( ntk, sim );
  return tts[0] == target;
}

} /* namespace */

int main( int argc, char** argv )
{
  uint32_t seed = 0xC0FFEEu;
  uint32_t reps = 200u;
  std::string aig_path;
  bool verbose = false;
  for ( int i = 1; i < argc; ++i )
  {
    std::string const a = argv[i];
    auto const eq = a.find( '=' );
    std::string const key = eq == std::string::npos ? a : a.substr( 0, eq );
    std::string const val = eq == std::string::npos ? "" : a.substr( eq + 1 );
    if ( key == "--seed" ) seed = std::stoul( val );
    else if ( key == "--reps" ) reps = std::stoul( val );
    else if ( key == "--aig" ) aig_path = val;
    else if ( key == "--verbose" ) verbose = true;
    else { std::cerr << "unknown option " << a << "\n"; return 1; }
  }

  std::mt19937 rng( seed );

  /* Two configurations.  `lut_size = 4` with the complete 4-input NPN database is the
   * exact-and-fast pairing; `lut_size = 6` needs an engine that covers six variables, for
   * which SOP factoring is the cheap choice.  Both are exercised because the cascade's
   * block width is the parameter most likely to break the decompArray walk. */
  /* the DEFAULT database is `xag_incomplete` and misses ~2/3 of the 4-input NPN
   * classes, which shows up here as block failures.  Ask for the complete one. */
  using npn4_t = xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete>;
  npn4_t npn4;
  sop_factoring<xag_network> sop;

  bool all_ok = true;
  uint32_t total_calls = 0, total_emitted = 0, total_acd = 0, total_bad = 0;

  /* `use_fallback` is exercised only with SOP factoring: it is the one of the two that can
   * legally be handed a function wider than its block size (see acd.hpp). */
  struct config { const char* name; uint32_t lut_size; bool use_npn; bool fallback; };
  config const configs[] = { { "lut4/xag_npn", 4u, true, false },
                             { "lut6/sop_factoring", 6u, false, true } };

  for ( auto const& cfg : configs )
  {
    for ( uint32_t n = 4; n <= 11; ++n )
    {
      outcome o;
      uint32_t eval_ok = 0;   /* acd_iface::evaluate said feasible                */
      acd_resynthesis_stats st;
      acd_resynthesis_params ps;
      ps.lut_size = cfg.lut_size;
      ps.use_fallback = cfg.fallback;

      /* the functor holds a reference to the stats block, so build it per (cfg, n) */
      acd_resynthesis<xag_network, npn4_t> resyn_npn( npn4, ps, &st );
      acd_resynthesis<xag_network, sop_factoring<xag_network>> resyn_sop( sop, ps, &st );

      for ( uint32_t r = 0; r < reps; ++r )
      {
        kitty::dynamic_truth_table target( n );
        kitty::create_random( target, rng() );

        /* the census used `evaluate`; the operator has to use `decompose`, which also
         * runs `compute_decomposition`.  Measure the gap between the two rather than
         * assuming there is none. */
        if ( n > cfg.lut_size )
        {
          std::vector<uint64_t> probe_bits( target._bits.begin(), target._bits.end() );
          uint32_t cost = 0;
          if ( acd_iface::evaluate( probe_bits.data(), n, cfg.lut_size, &cost ) >= 0 )
            ++eval_ok;
        }

        acd_resynthesis_stats const before = st;
        bool emitted = false, by_acd = false;
        bool ok;
        if ( cfg.use_npn )
          ok = check_one<xag_network>( resyn_npn, target, false, &emitted, &by_acd, before, st );
        else
          ok = check_one<xag_network>( resyn_sop, target, false, &emitted, &by_acd, before, st );

        ++o.calls;
        if ( emitted ) ++o.emitted;
        if ( by_acd ) ++o.accepted;
        if ( !ok ) ++o.mismatch;
      }

      total_calls += o.calls;
      total_emitted += o.emitted;
      total_acd += o.accepted;
      total_bad += o.mismatch;
      if ( o.mismatch != 0 )
        all_ok = false;

      fmt::print( "{:<20} n={:<3} calls={:<5} eval-ok={:<5} via-ACD={:<5} emitted={:<5} "
                  "blocks={:<6} widest-block={:<3} block-fail={:<4} (first at {:<2} fanins) "
                  "MISMATCH={}\n",
                  cfg.name, n, o.calls, eval_ok, o.accepted, o.emitted, st.num_blocks,
                  st.max_block_fanins, st.num_block_failures, st.first_failed_fanins,
                  o.mismatch );
      if ( verbose )
        fmt::print( "    stats: accepted={} declined={} trivial={}\n",
                    st.num_accepted, st.num_declined, st.num_trivial );
    }
  }

  /* ------------- real cut functions from an AIG, 7..11 variables ------------- */
  if ( !aig_path.empty() )
  {
    aig_network aig;
    if ( lorina::read_aiger( aig_path, aiger_reader( aig ) ) != lorina::return_code::success )
    {
      fmt::print( "FAIL: could not read {}\n", aig_path );
      return 4;
    }

    cut_enumeration_params cps;
    cps.cut_size = 11u;
    cps.cut_limit = 25u;
    auto const cuts = cut_enumeration<aig_network, true, check_wide_cut>( aig, cps );

    /* one collection per arity so the report shows where the coverage actually is */
    std::vector<kitty::dynamic_truth_table> pool;
    aig.foreach_gate( [&]( auto const& n ) {
      for ( auto const* cut : cuts.cuts( aig.node_to_index( n ) ) )
      {
        if ( cut->size() < 7u )
          continue;
        auto tt = cuts.truth_table( *cut );
        kitty::min_base_inplace( tt );
        uint32_t support = 0;
        for ( uint32_t v = 0; v < tt.num_vars(); ++v )
          if ( kitty::has_var( tt, v ) )
            ++support;
        if ( support < 7u || support > 11u )
          continue;
        pool.push_back( kitty::shrink_to( tt, support ) );
      }
    } );

    for ( uint32_t n = 7; n <= 11; ++n )
    {
      outcome o;
      acd_resynthesis_stats st;
      acd_resynthesis_params ps;
      ps.lut_size = 6u;
      ps.use_fallback = false;
      acd_resynthesis<xag_network, sop_factoring<xag_network>> resyn( sop, ps, &st );

      uint32_t taken = 0, eval_ok = 0;
      for ( auto const& tt : pool )
      {
        if ( tt.num_vars() != n || taken >= reps * 5u )
          continue;
        ++taken;
        {
          std::vector<uint64_t> probe_bits( tt._bits.begin(), tt._bits.end() );
          uint32_t cost = 0;
          if ( acd_iface::evaluate( probe_bits.data(), n, 6u, &cost ) >= 0 )
            ++eval_ok;
        }
        acd_resynthesis_stats const before = st;
        bool emitted = false, by_acd = false;
        bool const ok = check_one<xag_network>( resyn, tt, false, &emitted, &by_acd, before, st );
        ++o.calls;
        if ( emitted ) ++o.emitted;
        if ( by_acd ) ++o.accepted;
        if ( !ok ) ++o.mismatch;
      }

      total_calls += o.calls;
      total_emitted += o.emitted;
      total_acd += o.accepted;
      total_bad += o.mismatch;
      if ( o.mismatch != 0 )
        all_ok = false;

      fmt::print( "{:<20} n={:<3} calls={:<5} eval-ok={:<5} via-ACD={:<5} emitted={:<5} "
                  "blocks={:<6} widest-block={:<3} block-fail={:<4} (first at {:<2} fanins) "
                  "MISMATCH={}\n",
                  "real-cuts/lut6", n, o.calls, eval_ok, o.accepted, o.emitted, st.num_blocks,
                  st.max_block_fanins, st.num_block_failures, st.first_failed_fanins,
                  o.mismatch );
    }
  }

  /* ---------------- the sabotage control, which must fire ---------------- */
  uint32_t sabotage_cases = 0, sabotage_caught = 0;
  {
    acd_resynthesis_stats st;
    acd_resynthesis_params ps;
    ps.lut_size = 4u;
    acd_resynthesis<xag_network, npn4_t> resyn( npn4, ps, &st );

    for ( uint32_t n = 5; n <= 9; ++n )
    {
      for ( uint32_t r = 0; r < 40; ++r )
      {
        kitty::dynamic_truth_table target( n );
        kitty::create_random( target, rng() );
        acd_resynthesis_stats const before = st;
        bool emitted = false, by_acd = false;
        bool const ok = check_one<xag_network>( resyn, target, /*sabotage=*/true,
                                                &emitted, &by_acd, before, st );
        if ( !emitted )
          continue; /* nothing was emitted, so there was nothing to sabotage */
        ++sabotage_cases;
        if ( !ok )
          ++sabotage_caught;
      }
    }
  }

  fmt::print( "\ntotals: calls={} emitted={} via-ACD={} mismatches={}\n",
              total_calls, total_emitted, total_acd, total_bad );
  fmt::print( "sabotage control: {} of {} deliberately broken replacements were caught\n",
              sabotage_caught, sabotage_cases );

  bool const sabotage_ok = sabotage_cases > 0 && sabotage_caught == sabotage_cases;
  if ( !sabotage_ok )
  {
    fmt::print( "FAIL: the sabotage control did not fire -- the checker is vacuous and the "
                "{} clean results above mean nothing.\n", total_calls );
    return 2;
  }
  if ( total_acd == 0 )
  {
    fmt::print( "FAIL: ACD never produced a cascade, so only the fallback path was tested.\n" );
    return 3;
  }
  if ( !all_ok )
  {
    fmt::print( "FAIL: {} emitted replacements are not equivalent to their target.\n", total_bad );
    return 1;
  }
  fmt::print( "PASS: every emitted replacement is exhaustively equivalent, and the "
              "sabotage control fired on all {} cases.\n", sabotage_cases );
  return 0;
}
