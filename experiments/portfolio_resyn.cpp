/* portfolio_resyn: a PER-NODE PORTFOLIO resynthesis engine.
 *
 * Phase 2 (`experiments/2026-08-26-p2-acd-resyn`) established that no single resynthesis
 * engine dominates on the operator
 *
 *     AIG -> lut_map(K=8, area) -> 8-LUT network -> node_resynthesis<aig>(engine) -> AIG
 *
 * and that a per-BENCHMARK best-of-all oracle is worth -24.79 % AIG nodes / -7.05 % mapped
 * 6-LUTs against DSD alone -- but that selecting the arm by AIG size captures only -1.59 %
 * of it, because the AIG-best engine is the mapped-best on just 10 of 30 benchmarks.
 * Selection is the whole problem.
 *
 * This driver moves the selection inside the operator.  For EVERY node of the 8-LUT network
 * it runs all six candidate engines into six throwaway scratch AIGs, measures each
 * candidate three ways, picks one by a stated criterion, and replays the winner into the
 * destination AIG.  Every criterion is computable at resynthesis time, so all three modes
 * are DEPLOYABLE POLICIES, not oracles.
 *
 * Usage:
 *   portfolio_resyn <input> <output> [options]
 *
 *   <input>    .aig (binary AIGER) or .blif (k-LUT netlist, decomposed back to an AIG)
 *   <output>   .aig (default) or .blif, per --emit
 *
 * Options:
 *   --select=S     p-area | p-depth | p-lut                     (portfolio arms)
 *                  dsd | acd | shannon | sopf | bidec | direct  (single-engine baselines,
 *                  run through the IDENTICAL scratch-and-replay machinery so that the
 *                  baseline's runtime is the honest 1x reference and its instrumentation
 *                  is the same as a portfolio arm's)
 *   --map-k=K      LUT size of the intermediate k-LUT network (default 8, Phase 2's frozen
 *                  operating point)
 *   --cut-limit=N  priority-cut list size (default 12)
 *   --acd-lut-size=N   block size of the ACD cascade (default 6, Phase 2's operating point)
 *   --lut-cost-k=K     LUT size used for the p-lut local ruler (default 6, the endpoint's)
 *   --emit=aig|blif    what to write (default aig)
 *   --emit-k=K         LUT size when --emit=blif (default 6)
 *   --csv=PATH         dump one row per node_resynthesis call: every candidate's AIG size,
 *                      AIG depth and local LUT count, and which engine won
 *   --check            exhaustively simulate EVERY candidate (2^num_vars) against the
 *                      function it was asked to implement; any mismatch is fatal
 *   --sabotage-rate=P  (with --check) deliberately BREAK 1-in-P candidates, so that the
 *                      checker has to catch them.  The run exits non-zero if any injected
 *                      sabotage is missed, and also if none was injected -- a control that
 *                      never fires proves nothing.  A sabotage run WRITES NO OUTPUT: it
 *                      really does build wrong circuits, and must never be mistaken for a
 *                      scored run.
 *   --sabotage-mode=M  minterm (default) | invert
 *                      `minterm` resynthesises the candidate from a copy of the target with
 *                      ONE truth-table bit flipped, so the emitted block is wrong on exactly
 *                      one of the 2^n input patterns.  That is the weakest possible error and
 *                      is precisely what a sampling checker, or one that compares against the
 *                      wrong reference, would miss -- so it tests that the check is genuinely
 *                      exhaustive, not merely that it compares.
 *                      `invert` complements the emitted signal inside the scratch network
 *                      before it is measured, simulated and compared: a coarse break that
 *                      must also be caught.
 *   --sabotage-seed=S  seed for the sabotage RNG (default 1)
 *   --verbose
 *
 * THE SIX CANDIDATES, in the fixed order that also defines the tie-break priority:
 *
 *   0 dsd      dsd_resynthesis over the shared fallback
 *   1 acd      acd_resynthesis (lut_size 6, use_fallback) over the shared fallback
 *   2 shannon  the shared fallback alone
 *   3 sopf     sop_factoring
 *   4 bidec    bidecomposition_resynthesis
 *   5 direct   dsd_resynthesis over Shannon(4) + the aig_complete NPN database -- exactly
 *              what `convert_klut_to_graph` does per node, which is Phase 2's `direct` arm
 *
 * The shared fallback is `shannon_resynthesis` down to 4 variables then a COMPLETE 4-input
 * NPN database (xag_complete), identical to Phase 2's.  `direct` differs from `dsd` only in
 * that database (aig_complete), which is what makes it Phase 2's `direct`.
 *
 * THE THREE LOCAL RULERS, all measured on the candidate's own scratch block after
 * `cleanup_dangling` (which strashes it exactly as the replay into the destination will):
 *
 *   size   the AIG nodes the candidate emits for this node
 *   depth  the AIG depth of the emitted block, measured from its leaves
 *   luts   the 6-LUT count of an area-oriented `lut_map` of the emitted block, in isolation
 *
 * `luts` is the derived, ENGINE-UNIFORM equivalent of the cost ACD reports natively (its
 * cascade block count, also recorded in the CSV as `acd_blocks` for comparison).  A
 * per-engine structural LUT count is not definable for `sopf` or `bidec` at all, so a
 * uniform ruler is the only way the arms can be compared; mapping the block in isolation is
 * the closest computable thing in the endpoint's own units.  Note the consequence, which is
 * stated in the pre-registration rather than discovered afterwards: a block whose function
 * has <= 6 variables is one LUT for EVERY engine, so `p-lut` can only differ from `p-area`
 * on nodes with more than `--lut-cost-k` inputs.  The CSV records `num_vars` per node so
 * that fraction is recoverable.
 *
 * TIE-BREAK, deterministic by construction: each arm ranks candidates by a 4-tuple whose
 * last component is the engine index above, so no two candidates ever compare equal.
 *
 *   p-area   (size,  depth, luts,  engine)
 *   p-depth  (depth, size,  luts,  engine)
 *   p-lut    (luts,  size,  depth, engine)
 *
 * BUILD NOTE (inherited from acd_resyn.cpp): the whole ACD include chain must precede
 * anything that transitively pulls in ABC's `abc_namespaces.h`, which redefines
 * `ABC_NAMESPACE_CXX_HEADER_START` to `pabc` and lands the vendored decomposer in the wrong
 * namespace.
 *
 * \author agentic-synthesis
 */

#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>

#include <kitty/bit_operations.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

/* MUST come first -- see the build note above. */
#include <mockturtle/algorithms/node_resynthesis/acd.hpp>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/bidecomposition.hpp>
#include <mockturtle/algorithms/node_resynthesis/dsd.hpp>
#include <mockturtle/algorithms/node_resynthesis/shannon.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>

using namespace mockturtle;

namespace
{

constexpr uint32_t NUM_ENGINES = 6u;
char const* const ENGINE_NAME[NUM_ENGINES] = { "dsd", "acd", "shannon", "sopf", "bidec", "direct" };

enum class criterion
{
  area,
  depth,
  lut
};

/* the shared fallback: Shannon down to 4 variables, then a COMPLETE 4-input NPN database.
 * The default database is `xag_incomplete` and misses roughly two thirds of the classes,
 * failing silently -- that already cost the Phase 1 correctness gate a day. */
using npn_xag_t = xag_npn_resynthesis<aig_network, xag_network, xag_npn_db_kind::xag_complete>;
using npn_aig_t = xag_npn_resynthesis<aig_network, xag_network, xag_npn_db_kind::aig_complete>;
using fb_xag_t = shannon_resynthesis<aig_network, npn_xag_t>;
using fb_aig_t = shannon_resynthesis<aig_network, npn_aig_t>;

bool ends_with( std::string const& s, std::string const& suffix )
{
  return s.size() >= suffix.size() && s.compare( s.size() - suffix.size(), suffix.size(), suffix ) == 0;
}

bool read_input( std::string const& path, aig_network& aig )
{
  if ( ends_with( path, ".aig" ) || ends_with( path, ".aag" ) )
    return lorina::read_aiger( path, aiger_reader( aig ) ) == lorina::return_code::success;

  if ( ends_with( path, ".blif" ) )
  {
    klut_network klut;
    names_view<klut_network> named{ klut };
    if ( lorina::read_blif( path, blif_reader( named ) ) != lorina::return_code::success )
      return false;
    aig = convert_klut_to_graph<aig_network>( named );
    return true;
  }

  std::cerr << "[portfolio_resyn] unrecognised input extension: " << path << "\n";
  return false;
}

double secs_since( std::chrono::steady_clock::time_point t )
{
  return std::chrono::duration<double>( std::chrono::steady_clock::now() - t ).count();
}

/* ------------------------------------------------------------------ the portfolio ---- */

struct portfolio_params
{
  criterion crit{ criterion::lut };
  /* when set, only this engine is evaluated: the single-engine baseline arms, run through
   * the identical machinery so the comparison is not confounded by a second code path */
  std::optional<uint32_t> single{};
  uint32_t acd_lut_size{ 6u };
  uint32_t lut_cost_k{ 6u };
  bool check{ false };
  uint32_t sabotage_rate{ 0u };      /* 1-in-N, 0 = off */
  uint32_t sabotage_seed{ 1u };
  bool sabotage_minterm{ true };     /* minterm mode (default) vs invert mode */
};

struct portfolio_stats
{
  uint32_t calls{ 0u };
  uint32_t calls_wide{ 0u }; /* num_vars > lut_cost_k: the only nodes where p-lut can differ
                                from p-area other than through the shared tie-break */
  std::array<uint32_t, NUM_ENGINES> wins{};
  std::array<uint32_t, NUM_ENGINES> declines{};
  uint64_t sel_lut_total{ 0u };
  uint64_t sel_size_total{ 0u };
  uint64_t sel_depth_total{ 0u };
  /* the same three totals for each engine run alone over the same nodes -- free, because
   * every candidate is built anyway, and it makes the portfolio's headroom visible */
  std::array<uint64_t, NUM_ENGINES> engine_lut_total{};
  std::array<uint64_t, NUM_ENGINES> engine_size_total{};
  uint32_t check_calls{ 0u };
  uint32_t check_mismatch{ 0u };
  uint32_t sabotage_injected{ 0u };
  uint32_t sabotage_caught{ 0u };
  uint32_t no_candidate{ 0u };
};

struct candidate
{
  bool ok{ false };
  uint32_t size{ 0u };
  uint32_t depth{ 0u };
  uint32_t luts{ 0u };
  uint32_t acd_blocks{ 0u };
  aig_network net; /* num_vars PIs, exactly one PO */
};

class portfolio_resynthesis
{
public:
  portfolio_resynthesis( portfolio_params const& ps, portfolio_stats* pst,
                         acd_resynthesis_stats* acd_pst, std::ostream* csv )
      : _ps( ps ), _pst( pst ), _csv( csv ), _rng( ps.sabotage_seed ),
        _fb_xag( 4u, &_npn_xag ), _fb_aig( 4u, &_npn_aig ),
        _acd_ps( make_acd_params( ps.acd_lut_size ) ),
        _dsd( _fb_xag ), _acd( _fb_xag, _acd_ps, acd_pst ), _direct( _fb_aig ),
        _acd_pst( acd_pst )
  {
  }

  template<typename LeavesIterator, typename Fn>
  void operator()( aig_network& ntk, kitty::dynamic_truth_table const& function,
                   LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    std::vector<signal<aig_network>> leaves( begin, end );
    uint32_t const nv = function.num_vars();

    std::array<candidate, NUM_ENGINES> cands;

    if ( _ps.single.has_value() )
    {
      build_candidate( *_ps.single, function, nv, cands[*_ps.single] );
    }
    else
    {
      for ( uint32_t e = 0u; e < NUM_ENGINES; ++e )
        build_candidate( e, function, nv, cands[e] );
    }

    /* rank */
    int32_t best = -1;
    std::array<uint32_t, 4> best_key{};
    for ( uint32_t e = 0u; e < NUM_ENGINES; ++e )
    {
      if ( !cands[e].ok )
        continue;
      auto const k = rank_key( e, cands[e] );
      if ( best < 0 || k < best_key )
      {
        best = static_cast<int32_t>( e );
        best_key = k;
      }
    }

    if ( _pst != nullptr )
    {
      ++_pst->calls;
      if ( nv > _ps.lut_cost_k )
        ++_pst->calls_wide;
      for ( uint32_t e = 0u; e < NUM_ENGINES; ++e )
      {
        if ( cands[e].ok )
        {
          _pst->engine_lut_total[e] += cands[e].luts;
          _pst->engine_size_total[e] += cands[e].size;
        }
        else if ( !_ps.single.has_value() || *_ps.single == e )
        {
          ++_pst->declines[e];
        }
      }
    }

    if ( best < 0 )
    {
      /* every engine declined.  This must not happen -- the shared fallback is total for
       * any arity -- so it is fatal rather than silently papered over. */
      if ( _pst != nullptr )
        ++_pst->no_candidate;
      throw std::runtime_error( fmt::format(
          "no engine produced a candidate for a {}-variable node", nv ) );
    }

    auto const w = static_cast<uint32_t>( best );
    if ( _pst != nullptr )
    {
      ++_pst->wins[w];
      _pst->sel_lut_total += cands[w].luts;
      _pst->sel_size_total += cands[w].size;
      _pst->sel_depth_total += cands[w].depth;
    }

    if ( _csv != nullptr )
      write_csv_row( nv, cands, w );

    auto const outs = cleanup_dangling( cands[w].net, ntk, leaves.begin(), leaves.end() );
    fn( outs[0] );
  }

private:
  static acd_resynthesis_params make_acd_params( uint32_t lut_size )
  {
    acd_resynthesis_params p;
    p.lut_size = lut_size;
    p.max_num_vars = 11u;
    p.use_fallback = true; /* node_resynthesis needs total coverage */
    return p;
  }

  std::array<uint32_t, 4> rank_key( uint32_t e, candidate const& c ) const
  {
    switch ( _ps.crit )
    {
    case criterion::area:
      return { c.size, c.depth, c.luts, e };
    case criterion::depth:
      return { c.depth, c.size, c.luts, e };
    default:
      return { c.luts, c.size, c.depth, e };
    }
  }

  uint32_t count_luts( aig_network& net ) const
  {
    if ( net.num_gates() == 0u )
      return 0u;
    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = _ps.lut_cost_k;
    ps.area_oriented_mapping = true;
    klut_network m = lut_map<aig_network, true>( net, ps );
    return m.num_gates();
  }

  void build_candidate( uint32_t e, kitty::dynamic_truth_table const& f, uint32_t nv,
                        candidate& c ) const
  {
    /* the sabotage decision is taken HERE, before anything is built, so that a sabotaged
     * candidate is a genuinely wrong circuit all the way through measurement, simulation
     * and comparison -- not a value the comparator is handed after the fact. */
    bool sabotaged = false;
    kitty::dynamic_truth_table target = f;
    if ( _ps.check && _ps.sabotage_rate != 0u && nv > 0u &&
         ( _rng() % _ps.sabotage_rate ) == 0u )
    {
      sabotaged = true;
      if ( _ps.sabotage_minterm )
        kitty::flip_bit( target, _rng() % ( uint64_t{ 1 } << nv ) );
    }

    aig_network s;
    std::vector<signal<aig_network>> pis;
    pis.reserve( nv );
    for ( uint32_t i = 0u; i < nv; ++i )
      pis.push_back( s.create_pi() );

    bool got = false;
    signal<aig_network> out = s.get_constant( false );
    auto const cb = [&]( signal<aig_network> const& x ) {
      if ( !got )
      {
        out = x;
        got = true;
      }
      return false; /* one candidate per engine is enough; the portfolio ranks engines */
    };

    uint32_t const acd_blocks_before = _acd_pst != nullptr ? _acd_pst->num_blocks : 0u;
    uint32_t const acd_acc_before = _acd_pst != nullptr ? _acd_pst->num_accepted : 0u;

    try
    {
      switch ( e )
      {
      case 0u: _dsd( s, target, pis.begin(), pis.end(), cb ); break;
      case 1u: _acd( s, target, pis.begin(), pis.end(), cb ); break;
      case 2u: _fb_xag( s, target, pis.begin(), pis.end(), cb ); break;
      case 3u: _sopf( s, target, pis.begin(), pis.end(), cb ); break;
      case 4u: _bidec( s, target, pis.begin(), pis.end(), cb ); break;
      default: _direct( s, target, pis.begin(), pis.end(), cb ); break;
      }
    }
    catch ( std::exception const& )
    {
      got = false;
    }

    if ( !got )
      return;

    if ( e == 1u && _acd_pst != nullptr && _acd_pst->num_accepted != acd_acc_before )
      c.acd_blocks = _acd_pst->num_blocks - acd_blocks_before;

    if ( sabotaged && !_ps.sabotage_minterm )
      out = s.create_not( out );

    s.create_po( out );
    c.net = cleanup_dangling( s );
    c.size = c.net.num_gates();
    {
      depth_view<aig_network> d{ c.net };
      c.depth = d.depth();
    }
    c.luts = count_luts( c.net );
    c.ok = true;

    if ( _ps.check )
      run_check( e, f, nv, c, sabotaged );
  }

  /* Exhaustive 2^nv simulation of the emitted block against the function it was ASKED for
   * (`f`, always -- never the sabotaged target), plus the control that must fire.  Note
   * what makes this non-vacuous: a sabotaged candidate is a wrong circuit built from a
   * wrong specification, so `tts[0]` is genuinely computed from a different network.  A
   * checker that compared the block against whatever it was built from, or that sampled
   * instead of enumerating, would pass it. */
  void run_check( uint32_t e, kitty::dynamic_truth_table const& f, uint32_t nv,
                  candidate const& c, bool sabotaged ) const
  {
    if ( nv == 0u )
      return;

    default_simulator<kitty::dynamic_truth_table> sim( static_cast<unsigned>( nv ) );
    auto const tts = simulate<kitty::dynamic_truth_table>( c.net, sim );
    bool const match = ( tts[0] == f );

    if ( _pst != nullptr )
    {
      ++_pst->check_calls;
      if ( sabotaged )
        ++_pst->sabotage_injected;
    }

    if ( sabotaged )
    {
      if ( match )
        throw std::runtime_error( fmt::format(
            "SABOTAGE NOT CAUGHT: engine {} on a {}-variable node -- THE CHECKER IS BLIND",
            ENGINE_NAME[e], nv ) );
      if ( _pst != nullptr )
        ++_pst->sabotage_caught;
      return;
    }

    if ( !match )
    {
      if ( _pst != nullptr )
        ++_pst->check_mismatch;
      throw std::runtime_error( fmt::format(
          "CHECK FAILED: engine {} emitted a wrong function on a {}-variable node",
          ENGINE_NAME[e], nv ) );
    }
  }

  void write_csv_row( uint32_t nv, std::array<candidate, NUM_ENGINES> const& cands,
                      uint32_t w ) const
  {
    /* `call` is the ordinal of this node_resynthesis call in topological order.  The
     * 8-LUT network and its topo order are identical across arms, so the ordinal joins
     * rows across arms and is the per-node key for any later re-analysis. */
    std::string row = fmt::format( "{},{}", _csv_row++, nv );
    for ( uint32_t e = 0u; e < NUM_ENGINES; ++e )
    {
      if ( cands[e].ok )
        row += fmt::format( ",{},{},{}", cands[e].size, cands[e].depth, cands[e].luts );
      else
        row += ",-1,-1,-1";
    }
    row += fmt::format( ",{},{}\n", cands[1].acd_blocks, ENGINE_NAME[w] );
    *_csv << row;
  }

  portfolio_params _ps;
  portfolio_stats* _pst;
  std::ostream* _csv;
  mutable std::mt19937 _rng;
  mutable uint64_t _csv_row{ 0u };

  mutable npn_xag_t _npn_xag;
  mutable npn_aig_t _npn_aig;
  mutable fb_xag_t _fb_xag;
  mutable fb_aig_t _fb_aig;
  acd_resynthesis_params _acd_ps;
  mutable dsd_resynthesis<aig_network, fb_xag_t> _dsd;
  mutable acd_resynthesis<aig_network, fb_xag_t> _acd;
  mutable dsd_resynthesis<aig_network, fb_aig_t> _direct;
  mutable sop_factoring<aig_network> _sopf;
  mutable bidecomposition_resynthesis<aig_network> _bidec;
  acd_resynthesis_stats* _acd_pst;
};

} /* namespace */

int main( int argc, char** argv )
{
  if ( argc < 3 )
  {
    std::cerr << "usage: portfolio_resyn <input> <output> [options]\n";
    return 1;
  }

  std::string const input = argv[1];
  std::string const output = argv[2];

  std::string select = "p-lut";
  std::string emit = "aig";
  std::string csv_path;
  uint32_t map_k = 8u;
  uint32_t cut_limit = 12u;
  uint32_t emit_k = 6u;
  portfolio_params pps;
  bool verbose = false;

  for ( int i = 3; i < argc; ++i )
  {
    std::string const a = argv[i];
    auto const eq = a.find( '=' );
    std::string const key = eq == std::string::npos ? a : a.substr( 0, eq );
    std::string const val = eq == std::string::npos ? "" : a.substr( eq + 1 );

    if ( key == "--select" ) select = val;
    else if ( key == "--emit" ) emit = val;
    else if ( key == "--csv" ) csv_path = val;
    else if ( key == "--map-k" ) map_k = std::stoul( val );
    else if ( key == "--cut-limit" ) cut_limit = std::stoul( val );
    else if ( key == "--emit-k" ) emit_k = std::stoul( val );
    else if ( key == "--acd-lut-size" ) pps.acd_lut_size = std::stoul( val );
    else if ( key == "--lut-cost-k" ) pps.lut_cost_k = std::stoul( val );
    else if ( key == "--check" ) pps.check = true;
    else if ( key == "--sabotage-rate" ) pps.sabotage_rate = std::stoul( val );
    else if ( key == "--sabotage-seed" ) pps.sabotage_seed = std::stoul( val );
    else if ( key == "--sabotage-mode" )
    {
      if ( val == "minterm" ) pps.sabotage_minterm = true;
      else if ( val == "invert" ) pps.sabotage_minterm = false;
      else { std::cerr << "[portfolio_resyn] unknown --sabotage-mode: " << val << "\n"; return 1; }
    }
    else if ( key == "--verbose" ) verbose = true;
    else
    {
      std::cerr << "[portfolio_resyn] unknown option: " << a << "\n";
      return 1;
    }
  }

  if ( select == "p-area" ) pps.crit = criterion::area;
  else if ( select == "p-depth" ) pps.crit = criterion::depth;
  else if ( select == "p-lut" ) pps.crit = criterion::lut;
  else
  {
    bool found = false;
    for ( uint32_t e = 0u; e < NUM_ENGINES; ++e )
    {
      if ( select == ENGINE_NAME[e] )
      {
        pps.single = e;
        found = true;
        break;
      }
    }
    if ( !found )
    {
      std::cerr << "[portfolio_resyn] unknown --select: " << select << "\n";
      return 1;
    }
  }

  if ( pps.acd_lut_size < 2u || pps.acd_lut_size > 6u )
  {
    std::cerr << "[portfolio_resyn] --acd-lut-size must be in [2,6]\n";
    return 1;
  }
  if ( emit != "aig" && emit != "blif" )
  {
    std::cerr << "[portfolio_resyn] unknown emit mode: " << emit << "\n";
    return 1;
  }
  if ( pps.sabotage_rate != 0u && !pps.check )
  {
    std::cerr << "[portfolio_resyn] --sabotage-rate requires --check\n";
    return 1;
  }

  aig_network aig_in;
  if ( !read_input( input, aig_in ) )
  {
    std::cerr << "[portfolio_resyn] could not read " << input << "\n";
    return 1;
  }
  aig_in = cleanup_dangling( aig_in );

  uint32_t const pi_in = aig_in.num_pis();
  uint32_t const po_in = aig_in.num_pos();
  uint32_t const aig_in_size = aig_in.num_gates();
  uint32_t aig_in_depth = 0;
  {
    depth_view<aig_network> d{ aig_in };
    aig_in_depth = d.depth();
  }

  auto const t_start = std::chrono::steady_clock::now();

  /* ---- the shared intermediate k-LUT network ------------------------------------- */
  auto const t_map = std::chrono::steady_clock::now();
  lut_map_params mps;
  mps.cut_enumeration_ps.cut_size = map_k;
  mps.cut_enumeration_ps.cut_limit = cut_limit;
  mps.area_oriented_mapping = true;
  mps.verbose = verbose;
  klut_network klut = lut_map<aig_network, true>( aig_in, mps );
  double const secs_map = secs_since( t_map );

  uint32_t const klut_size = klut.num_gates();
  uint32_t klut_max_fanin = 0u;
  klut.foreach_gate( [&]( auto const& n ) {
    klut_max_fanin = std::max( klut_max_fanin, klut.fanin_size( n ) );
  } );

  /* ---- resynthesis ---------------------------------------------------------------- */
  std::ofstream csv;
  if ( !csv_path.empty() )
  {
    csv.open( csv_path );
    if ( !csv )
    {
      std::cerr << "[portfolio_resyn] could not open " << csv_path << " for writing\n";
      return 1;
    }
    csv << "call,num_vars";
    for ( uint32_t e = 0u; e < NUM_ENGINES; ++e )
      csv << fmt::format( ",{0}_size,{0}_depth,{0}_luts", ENGINE_NAME[e] );
    csv << ",acd_blocks,winner\n";
  }

  portfolio_stats pst;
  acd_resynthesis_stats acd_st;
  aig_network aig_out;

  double secs_resyn = 0.0;
  try
  {
    portfolio_resynthesis resyn( pps, &pst, &acd_st, csv_path.empty() ? nullptr : &csv );
    auto const t_resyn = std::chrono::steady_clock::now();
    aig_out = node_resynthesis<aig_network>( klut, resyn );
    secs_resyn = secs_since( t_resyn );
  }
  catch ( std::exception const& e )
  {
    std::cerr << "[portfolio_resyn] FATAL: " << e.what() << "\n";
    return 3;
  }

  if ( csv.is_open() )
    csv.close();

  aig_out = cleanup_dangling( aig_out );

  if ( acd_st.num_block_failures != 0u )
  {
    std::cerr << fmt::format( "[portfolio_resyn] FATAL: {} ACD blocks could not be realised "
                              "(first at {} fanins)\n",
                              acd_st.num_block_failures, acd_st.first_failed_fanins );
    return 1;
  }

  if ( aig_out.num_pis() != pi_in || aig_out.num_pos() != po_in )
  {
    std::cerr << fmt::format( "[portfolio_resyn] FATAL: interface changed, PI {}->{} PO {}->{}\n",
                              pi_in, aig_out.num_pis(), po_in, aig_out.num_pos() );
    return 1;
  }

  uint32_t const aig_out_size = aig_out.num_gates();
  uint32_t aig_out_depth = 0;
  {
    depth_view<aig_network> d{ aig_out };
    aig_out_depth = d.depth();
  }

  /* A sabotage run really does build wrong circuits.  It writes nothing, so it can never
   * be scored by accident and can never be mistaken for a result. */
  if ( pps.sabotage_rate != 0u )
  {
    std::cerr << "[portfolio_resyn] sabotage control run: no output written\n";
  }
  else if ( emit == "aig" )
  {
    write_aiger( aig_out, output );
  }
  else
  {
    lut_map_params eps;
    eps.cut_enumeration_ps.cut_size = emit_k;
    eps.area_oriented_mapping = true;
    klut_network mapped = lut_map<aig_network, true>( aig_out, eps );

    uint32_t widest = 0;
    mapped.foreach_gate( [&]( auto const& n ) {
      widest = std::max( widest, mapped.fanin_size( n ) );
    } );
    if ( widest > emit_k )
    {
      std::cerr << fmt::format( "[portfolio_resyn] FATAL: widest LUT has {} inputs, limit {}\n",
                                widest, emit_k );
      return 1;
    }
    write_blif( mapped, output );
  }

  double const secs_total = secs_since( t_start );

  std::string wins, declines, eng_luts, eng_sizes;
  for ( uint32_t e = 0u; e < NUM_ENGINES; ++e )
  {
    wins += fmt::format( " win_{}={}", ENGINE_NAME[e], pst.wins[e] );
    declines += fmt::format( " dec_{}={}", ENGINE_NAME[e], pst.declines[e] );
    eng_luts += fmt::format( " lutsum_{}={}", ENGINE_NAME[e], pst.engine_lut_total[e] );
    eng_sizes += fmt::format( " sizesum_{}={}", ENGINE_NAME[e], pst.engine_size_total[e] );
  }

  fmt::print( stderr,
              "[portfolio_resyn] select={} map_k={} lut_cost_k={} pi={} po={} "
              "aig_in_size={} aig_in_depth={} klut_size={} klut_max_fanin={} "
              "aig_out_size={} aig_out_depth={} "
              "calls={} calls_wide={} no_candidate={}"
              "{}{}"
              " sel_lut_total={} sel_size_total={} sel_depth_total={}"
              "{}{}"
              " acd_accepted={} acd_declined={} acd_trivial={} acd_blocks={} "
              "acd_block_failures={} acd_max_block_fanins={} "
              "check_calls={} check_mismatch={} sabotage_injected={} sabotage_caught={} "
              "secs_map={:.3f} secs_resyn={:.3f} secs_total={:.3f}\n",
              select, map_k, pps.lut_cost_k, pi_in, po_in,
              aig_in_size, aig_in_depth, klut_size, klut_max_fanin,
              aig_out_size, aig_out_depth,
              pst.calls, pst.calls_wide, pst.no_candidate,
              wins, declines,
              pst.sel_lut_total, pst.sel_size_total, pst.sel_depth_total,
              eng_luts, eng_sizes,
              acd_st.num_accepted, acd_st.num_declined, acd_st.num_trivial, acd_st.num_blocks,
              acd_st.num_block_failures, acd_st.max_block_fanins,
              pst.check_calls, pst.check_mismatch, pst.sabotage_injected, pst.sabotage_caught,
              secs_map, secs_resyn, secs_total );

  if ( pst.check_mismatch != 0u )
  {
    std::cerr << "[portfolio_resyn] FATAL: correctness check reported mismatches\n";
    return 3;
  }
  if ( pps.sabotage_rate != 0u )
  {
    if ( pst.sabotage_injected == 0u )
    {
      std::cerr << "[portfolio_resyn] FATAL: sabotage control injected nothing -- "
                   "the control is vacuous\n";
      return 4;
    }
    if ( pst.sabotage_caught != pst.sabotage_injected )
    {
      std::cerr << fmt::format( "[portfolio_resyn] FATAL: sabotage control missed {} of {}\n",
                                pst.sabotage_injected - pst.sabotage_caught,
                                pst.sabotage_injected );
      return 4;
    }
    /* a sabotage run is a control, not a result: never let it be mistaken for one */
    std::cerr << fmt::format( "[portfolio_resyn] SABOTAGE CONTROL PASSED: {}/{} caught\n",
                              pst.sabotage_caught, pst.sabotage_injected );
  }

  return 0;
}
