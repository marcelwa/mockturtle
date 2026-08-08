/* mt_flow: a harness-facing driver for the agentic-synthesis project.
 *
 * Usage:
 *   mt_flow <input> <output.blif> [options]
 *
 *   <input>       .aig (binary AIGER), .v (structural Verilog) or .blif (k-LUT netlist,
 *                 e.g. one of our champions -- it is decomposed back into an AIG).
 *   <output.blif> a 6-LUT mapped BLIF, ready for the harness legality gate.
 *
 * Options:
 *   --flow=<a,b,c>     comma-separated optimisation ops, applied left to right.
 *   --rounds=N         repeat the whole flow N times (default 1); stops early on no gain.
 *   --mig-flow=<...>   ops used inside the `mig` excursion op.
 *   --k=N              LUT size for the final mapping (default 6).
 *   --cut-limit=N      cut limit for the final mapping (default 8).
 *   --map=area|delay|sop|esop|mffc   final mapping style (default area).
 *   --relax=N          required-delay relaxation in % for the final mapping.
 *   --max-pis=N        resubstitution window inputs (default 8).
 *   --max-inserts=N    resubstitution insertion limit (default 2).
 *   --max-divisors=N   resubstitution divisor limit (default 150).
 *   --seed=N           random seed where applicable.
 *   --verbose          progress on stderr.
 *
 * SPFD support selection (op `rspfd`, see resyn_engines/spfd_resyn.hpp):
 *   --spfd-k=N         max support size K explored by the covering process (default 7).
 *   --spfd-samples=N   supports sampled per resynthesis call, S (default 10). 0 = off.
 *   --spfd-beta=F      inverse temperature on the normalised remaining-edge count;
 *                      negative means pure greedy support selection (default 5).
 *   --spfd-max-divs=N  divisors entering the covering process (default 150).
 *   --spfd-onfail      only sample supports when the plain engine found nothing.
 *
 * XAG tail (the information-graph line of work is an XAIG method):
 *   --xag-flow=<ops>   after the AIG flow, convert to an XAG, run these ops, and map to
 *                      LUTs from the XAG. Ops: xrw (NPN-4 XAG rewrite), xrs
 *                      (sim_resubstitution), xspfd (SPFD sim_resubstitution), xb (SOP
 *                      rebalancing). `aig_out=` in the summary line then reports XAG nodes.
 *
 * Design-space exploration (explorer.hpp / deepsyn -- basin hopping):
 *   --explore=<engine>   run an explorer before `--flow`. Engines:
 *                          aig     deepsyn_aig      (needs ENABLE_ABC)
 *                          migv1   deepsyn_mig_v1   (needs ENABLE_ABC)
 *                          migv2   deepsyn_mig_v2   (needs ENABLE_ABC)
 *                          migd    deepsyn_mig_depth(needs ENABLE_ABC)
 *                          mig     explore_mig      (ABC-free)
 *   --explore-cost=size|lut   cost the explorer minimises. `size` is the shipped
 *                          behaviour (AND/MAJ gate count); `lut` is the *harness metric*
 *                          -- the number of k-LUTs the final mapper would emit -- which
 *                          makes the search optimise what we are actually scored on.
 *   --restarts=N         explorer_params::num_restarts (default 1).
 *   --explore-timeout=S  seconds per restart (default 300).
 *   --steps=N            max steps per restart (default 100000).
 *   --steps-no-impr=N    give up a restart after N steps without improvement.
 *   --compress-per-step=N  compressing scripts per step (default 3).
 *
 * AIG ops:
 *   b     aig_balance (level-minimising)
 *   bf    aig_balance (fast, no level minimisation)
 *   sopb  SOP rebalancing (depth-oriented, changes structure a lot)
 *   rw    rewrite, AIG NPN-4 database
 *   rwz   rewrite, zero-gain moves allowed (non-monotone)
 *   rwd   rewrite with don't cares
 *   rf    refactoring via SOP factoring
 *   rfz   refactoring, zero-gain moves allowed
 *   rs    aig_resubstitution
 *   rs2   aig_resubstitution2 (with the resub engine of aig_resub.hpp)
 *   rsim  sim_resubstitution (simulation-guided + SAT validation)
 *   rspfd sim_resubstitution with SPFD/information-graph statistical support selection
 *         (Costamagna et al.). Tuned by --spfd-*; --spfd-samples=0 makes it identical to
 *         `rsim` and is the matched control arm.
 *   wr    window_rewriting
 *   wrd   window_rewriting with don't cares
 *   fr    functional_reduction (SAT-based structural hashing across the network)
 *   mig   MIG excursion: AIG -> MIG, run --mig-flow, MIG -> AIG
 *
 * MIG ops (inside --mig-flow):
 *   mrs   mig_resubstitution
 *   mrs2  mig_resubstitution2
 *   mad   mig_algebraic_depth_rewriting
 *   mrw   rewrite with the MIG NPN-4 database
 *   mmap  exact-library mapping onto MIG (area-oriented)
 *   msopb SOP rebalancing on the MIG
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>
#include <lorina/verilog.hpp>

#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/aig_resub.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/explorer.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/mig_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/spfd_resub.hpp>
#include <mockturtle/algorithms/window_rewriting.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/names_view.hpp>

using namespace mockturtle;

namespace
{

bool g_verbose = false;

struct options
{
  std::string flow = "b,rs,rw,rf,rs,rw,rs";
  std::string mig_flow = "mrw,mrs,mad";
  uint32_t rounds = 1u;
  uint32_t k = 6u;
  uint32_t cut_limit = 8u;
  std::string map_style = "area";
  uint32_t relax = 0u;
  uint32_t max_pis = 8u;
  uint32_t max_inserts = 2u;
  uint32_t max_divisors = 150u;
  uint32_t seed = 1u;
  std::string explore = "";
  std::string explore_cost = "size";
  uint32_t restarts = 1u;
  uint32_t explore_timeout = 300u;
  uint32_t steps = 100000u;
  uint32_t steps_no_impr = 1000000u;
  uint32_t compress_per_step = 3u;
  uint32_t spfd_k = 7u;
  uint32_t spfd_samples = 10u;
  double spfd_beta = 5.0;
  uint32_t spfd_max_divs = 150u;
  bool spfd_onfail = false;
  bool spfd_diagnose = false;
  std::string xag_flow = "";
};

std::vector<std::string> split( std::string const& s, char sep )
{
  std::vector<std::string> out;
  std::string cur;
  std::istringstream is( s );
  while ( std::getline( is, cur, sep ) )
  {
    /* trim */
    auto b = cur.find_first_not_of( " \t" );
    if ( b == std::string::npos )
      continue;
    auto e = cur.find_last_not_of( " \t" );
    out.push_back( cur.substr( b, e - b + 1 ) );
  }
  return out;
}

bool ends_with( std::string const& s, std::string const& suffix )
{
  return s.size() >= suffix.size() && s.compare( s.size() - suffix.size(), suffix.size(), suffix ) == 0;
}

void log( std::string const& msg )
{
  if ( g_verbose )
    std::cerr << "[mt_flow] " << msg << "\n";
}

/* ---------------------------------------------------------------- reading */

bool read_input( std::string const& path, aig_network& aig )
{
  if ( ends_with( path, ".aig" ) || ends_with( path, ".aag" ) )
  {
    return lorina::read_aiger( path, aiger_reader( aig ) ) == lorina::return_code::success;
  }
  if ( ends_with( path, ".v" ) || ends_with( path, ".verilog" ) )
  {
    return lorina::read_verilog( path, verilog_reader( aig ) ) == lorina::return_code::success;
  }
  if ( ends_with( path, ".blif" ) )
  {
    /* A mapped k-LUT netlist -- one of our own champions, say. Decompose it back into an
     * AIG so that the AIG-level algorithms have something to work on. This is what makes
     * `start_from: champion:size` usable from this driver. */
    klut_network klut;
    names_view<klut_network> named{ klut };
    if ( lorina::read_blif( path, blif_reader( named ) ) != lorina::return_code::success )
      return false;
    aig = convert_klut_to_graph<aig_network>( named );
    return true;
  }
  std::cerr << "[mt_flow] unrecognised input extension: " << path << "\n";
  return false;
}

/* ------------------------------------------------------------- MIG domain */

void run_mig_op( mig_network& mig, std::string const& op, options const& opts )
{
  if ( op == "mrs" || op == "mrs2" )
  {
    resubstitution_params ps;
    ps.max_pis = opts.max_pis;
    ps.max_inserts = opts.max_inserts;
    ps.max_divisors = opts.max_divisors;
    depth_view depth_mig{ mig };
    fanout_view fanout_mig{ depth_mig };
    if ( op == "mrs" )
      mig_resubstitution( fanout_mig, ps );
    else
      mig_resubstitution2( fanout_mig, ps );
    mig = cleanup_dangling( mig );
  }
  else if ( op == "mad" )
  {
    depth_view depth_mig{ mig };
    mig_algebraic_depth_rewriting( depth_mig );
    mig = cleanup_dangling( mig );
  }
  else if ( op == "mrw" || op == "mmap" )
  {
    mig_npn_resynthesis resyn{ true };
    exact_library_params eps;
    exact_library<mig_network> lib( resyn, eps );
    if ( op == "mrw" )
    {
      rewrite_params ps;
      rewrite( mig, lib, ps );
      mig = cleanup_dangling( mig );
    }
    else
    {
      map_params mps;
      mps.skip_delay_round = true;
      mps.required_time = std::numeric_limits<double>::max();
      mig = map( mig, lib, mps );
    }
  }
  else if ( op == "msopb" )
  {
    sop_rebalancing<mig_network> balance_fn;
    balancing_params bps;
    bps.cut_enumeration_ps.cut_size = 6u;
    mig = balancing( mig, { balance_fn }, bps );
  }
  else
  {
    std::cerr << "[mt_flow] unknown MIG op: " << op << "\n";
  }
}

/* ------------------------------------------------------------- AIG domain */

void run_aig_op( aig_network& aig, std::string const& op, options const& opts )
{
  auto const before = aig.num_gates();

  if ( op == "b" || op == "bf" )
  {
    aig_balancing_params ps;
    ps.minimize_levels = ( op == "b" );
    aig_balance( aig, ps );
  }
  else if ( op == "sopb" )
  {
    sop_rebalancing<aig_network> balance_fn;
    balancing_params bps;
    bps.cut_enumeration_ps.cut_size = opts.k;
    aig = balancing( aig, { balance_fn }, bps );
  }
  else if ( op == "rw" || op == "rwz" || op == "rwd" )
  {
    xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
    exact_library_params eps;
    eps.compute_dc_classes = ( op == "rwd" );
    exact_library<aig_network> lib( resyn, eps );
    rewrite_params ps;
    ps.allow_zero_gain = ( op == "rwz" );
    ps.use_dont_cares = ( op == "rwd" );
    rewrite( aig, lib, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "rf" || op == "rfz" )
  {
    sop_factoring<aig_network> resyn;
    refactoring_params ps;
    ps.max_pis = 10u;
    ps.allow_zero_gain = ( op == "rfz" );
    refactoring( aig, resyn, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "rs" || op == "rs2" )
  {
    resubstitution_params ps;
    ps.max_pis = opts.max_pis;
    ps.max_inserts = opts.max_inserts;
    ps.max_divisors = opts.max_divisors;
    depth_view depth_aig{ aig };
    fanout_view fanout_aig{ depth_aig };
    if ( op == "rs" )
      aig_resubstitution( fanout_aig, ps );
    else
      aig_resubstitution2( fanout_aig, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "rsim" )
  {
    resubstitution_params ps;
    ps.max_pis = opts.max_pis;
    ps.max_inserts = opts.max_inserts;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();
    ps.random_seed = opts.seed;
    sim_resubstitution( aig, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "rspfd" )
  {
    resubstitution_params ps;
    ps.max_pis = opts.max_pis;
    ps.max_inserts = opts.max_inserts;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();
    ps.random_seed = opts.seed;
    ps.verbose = g_verbose;
    auto& sp = spfd_global_params();
    sp.max_support = opts.spfd_k;
    sp.num_supports = opts.spfd_samples;
    sp.beta = opts.spfd_beta;
    sp.max_divisors = opts.spfd_max_divs;
    sp.seed = opts.seed;
    sp.only_on_fail = opts.spfd_onfail;
    sp.diagnose = opts.spfd_diagnose;
    spfd_sim_resubstitution( aig, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "wr" || op == "wrd" )
  {
    window_rewriting_params ps;
    ps.cut_size = 6u;
    ps.num_levels = 5u;
    ps.filter_cyclic_substitutions = true;
    ps.use_dont_cares = ( op == "wrd" );
    window_rewriting( aig, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "fr" )
  {
    functional_reduction_params ps;
    functional_reduction( aig, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "mig" )
  {
    mig_network mig = cleanup_dangling<aig_network, mig_network>( aig );
    for ( auto const& mop : split( opts.mig_flow, ',' ) )
    {
      auto const mbefore = mig.num_gates();
      run_mig_op( mig, mop, opts );
      log( fmt::format( "  mig/{}: {} -> {}", mop, mbefore, mig.num_gates() ) );
    }
    aig = cleanup_dangling<mig_network, aig_network>( mig );
  }
  else
  {
    std::cerr << "[mt_flow] unknown AIG op: " << op << "\n";
    return;
  }

  log( fmt::format( "{}: {} -> {} gates", op, before, aig.num_gates() ) );
}

/* ------------------------------------------------------------- XAG domain */

/* The information-graph line of work operates on XAIGs, not AIGs, and that turns out to
 * matter: the supports SPFD selection finds but unateness-based selection misses are
 * XOR-dominant, and on an AIG an XOR costs three nodes, so those supports are never
 * payable inside a resubstitution budget derived from a small MFFC. On an XAG they cost
 * one node. `--xag-flow` runs the tail of the flow, and the LUT mapping, on an XAG. */
void run_xag_op( xag_network& xag, std::string const& op, options const& opts )
{
  auto const before = xag.num_gates();

  if ( op == "xrw" )
  {
    xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete> resyn;
    exact_library_params eps;
    exact_library<xag_network> lib( resyn, eps );
    rewrite_params ps;
    rewrite( xag, lib, ps );
    xag = cleanup_dangling( xag );
  }
  else if ( op == "xrs" || op == "xspfd" )
  {
    resubstitution_params ps;
    ps.max_pis = opts.max_pis;
    ps.max_inserts = opts.max_inserts;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();
    ps.random_seed = opts.seed;
    ps.verbose = g_verbose;
    if ( op == "xrs" )
    {
      sim_resubstitution( xag, ps );
    }
    else
    {
      auto& sp = spfd_global_params();
      sp.max_support = opts.spfd_k;
      sp.num_supports = opts.spfd_samples;
      sp.beta = opts.spfd_beta;
      sp.max_divisors = opts.spfd_max_divs;
      sp.seed = opts.seed;
      sp.only_on_fail = opts.spfd_onfail;
    sp.diagnose = opts.spfd_diagnose;
      spfd_sim_resubstitution( xag, ps );
    }
    xag = cleanup_dangling( xag );
  }
  else if ( op == "xb" )
  {
    sop_rebalancing<xag_network> balance_fn;
    balancing_params bps;
    bps.cut_enumeration_ps.cut_size = opts.k;
    xag = balancing( xag, { balance_fn }, bps );
  }
  else
  {
    std::cerr << "[mt_flow] unknown XAG op: " << op << "\n";
    return;
  }

  log( fmt::format( "xag/{}: {} -> {} gates", op, before, xag.num_gates() ) );
}

/* -------------------------------------------------------------- mapping */

template<class Ntk>
klut_network map_to_luts( Ntk const& aig, options const& opts )
{
  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = opts.k;
  ps.cut_enumeration_ps.cut_limit = opts.cut_limit;
  ps.recompute_cuts = true;
  ps.cut_expansion = true;
  ps.relax_required = opts.relax;

  if ( opts.map_style == "area" )
  {
    ps.area_oriented_mapping = true;
  }
  else if ( opts.map_style == "delay" )
  {
    ps.area_oriented_mapping = false;
  }
  else if ( opts.map_style == "sop" )
  {
    ps.area_oriented_mapping = false;
    ps.sop_balancing = true;
  }
  else if ( opts.map_style == "esop" )
  {
    ps.area_oriented_mapping = false;
    ps.esop_balancing = true;
  }
  else if ( opts.map_style == "mffc" )
  {
    ps.area_oriented_mapping = true;
    ps.collapse_mffcs = true;
  }
  else
  {
    std::cerr << "[mt_flow] unknown --map style: " << opts.map_style << "\n";
  }

  return lut_map( aig, ps );
}

/* ------------------------------------------------------ explorer cost fns */

/* The harness counts a `.names` block as a LUT only when it has >= 2 inputs, and takes the
 * depth over exactly those blocks. Buffers and inverters emitted by the mapper are free.
 * Reproducing that rule here is what makes `--explore-cost=lut` optimise the scored metric
 * rather than a proxy for it. */
uint32_t count_scored_luts( klut_network const& klut )
{
  uint32_t n = 0u;
  klut.foreach_gate( [&]( auto const& g ) {
    if ( klut.fanin_size( g ) >= 2u )
      ++n;
  } );
  return n;
}

/* Number of k-LUTs the final area-oriented mapper would emit for `ntk`. Defined for any
 * network the explorer works on (AIG or MIG); `lut_map` needs a mutable network, so the
 * candidate is cloned. */
template<class Ntk>
uint32_t mapped_lut_cost( Ntk const& ntk, options const& opts )
{
  Ntk copy = ntk.clone();
  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = opts.k;
  ps.cut_enumeration_ps.cut_limit = opts.cut_limit;
  ps.recompute_cuts = true;
  ps.cut_expansion = true;
  ps.area_oriented_mapping = true;
  return count_scored_luts( lut_map( copy, ps ) );
}

template<class Ntk>
cost_fn_t<Ntk> explorer_cost( options const& opts )
{
  if ( opts.explore_cost == "lut" )
    return [&opts]( Ntk const& ntk ) { return mapped_lut_cost<Ntk>( ntk, opts ); };
  return size_cost_fn<Ntk>;
}

/* Run one of the explorer.hpp basin-hopping engines on `aig`. Returns false if the
 * requested engine is unavailable in this build. */
bool run_explorer( aig_network& aig, options const& opts )
{
  explorer_params eps;
  eps.num_restarts = opts.restarts;
  eps.random_seed = opts.seed;
  eps.max_steps = opts.steps;
  eps.max_steps_no_impr = opts.steps_no_impr;
  eps.compressing_scripts_per_step = opts.compress_per_step;
  eps.timeout = opts.explore_timeout;
  eps.verbose = g_verbose;

  if ( opts.explore == "mig" )
  {
    mig_network mig = cleanup_dangling<aig_network, mig_network>( aig );
    mig = explore_mig( mig, eps, explorer_cost<mig_network>( opts ) );
    aig = cleanup_dangling<mig_network, aig_network>( mig );
    return true;
  }
#ifdef ENABLE_ABC
  if ( opts.explore == "aig" )
  {
    aig = deepsyn_aig( aig, eps, explorer_cost<aig_network>( opts ) );
    return true;
  }
  if ( opts.explore == "migv1" || opts.explore == "migv2" || opts.explore == "migd" )
  {
    mig_network mig = cleanup_dangling<aig_network, mig_network>( aig );
    if ( opts.explore == "migv1" )
      mig = deepsyn_mig_v1( mig, eps, explorer_cost<mig_network>( opts ) );
    else if ( opts.explore == "migv2" )
      mig = deepsyn_mig_v2( mig, eps, explorer_cost<mig_network>( opts ) );
    else
      mig = deepsyn_mig_depth( mig, eps );
    aig = cleanup_dangling<mig_network, aig_network>( mig );
    return true;
  }
  std::cerr << "[mt_flow] unknown --explore engine: " << opts.explore << "\n";
#else
  std::cerr << "[mt_flow] --explore=" << opts.explore
            << " needs an ENABLE_ABC build (lib/abc_static/libabc.a); only "
               "--explore=mig is available here\n";
#endif
  return false;
}

} // namespace

int main( int argc, char** argv )
{
  if ( argc < 3 )
  {
    std::cerr << "usage: mt_flow <input.aig|.v|.blif> <output.blif> [--flow=...] "
                 "[--rounds=N] [--k=N] [--cut-limit=N] [--map=area|delay|sop|esop|mffc] "
                 "[--relax=N] [--max-pis=N] [--max-inserts=N] [--max-divisors=N] "
                 "[--mig-flow=...] [--seed=N] [--verbose] "
                 "[--explore=aig|mig|migv1|migv2|migd] [--explore-cost=size|lut] "
                 "[--restarts=N] [--explore-timeout=SEC] [--steps=N] [--steps-no-impr=N] "
                 "[--compress-per-step=N] [--spfd-k=N] [--spfd-samples=N] [--spfd-beta=F] "
                 "[--spfd-max-divs=N] [--spfd-onfail]\n";
    return 1;
  }

  std::string const input = argv[1];
  std::string const output = argv[2];
  options opts;

  for ( int i = 3; i < argc; ++i )
  {
    std::string a = argv[i];
    auto const eq = a.find( '=' );
    std::string key = eq == std::string::npos ? a : a.substr( 0, eq );
    std::string val = eq == std::string::npos ? "" : a.substr( eq + 1 );

    if ( key == "--flow" ) opts.flow = val;
    else if ( key == "--mig-flow" ) opts.mig_flow = val;
    else if ( key == "--rounds" ) opts.rounds = std::stoul( val );
    else if ( key == "--k" ) opts.k = std::stoul( val );
    else if ( key == "--cut-limit" ) opts.cut_limit = std::stoul( val );
    else if ( key == "--map" ) opts.map_style = val;
    else if ( key == "--relax" ) opts.relax = std::stoul( val );
    else if ( key == "--max-pis" ) opts.max_pis = std::stoul( val );
    else if ( key == "--max-inserts" ) opts.max_inserts = std::stoul( val );
    else if ( key == "--max-divisors" ) opts.max_divisors = std::stoul( val );
    else if ( key == "--seed" ) opts.seed = std::stoul( val );
    else if ( key == "--explore" ) opts.explore = val;
    else if ( key == "--explore-cost" ) opts.explore_cost = val;
    else if ( key == "--restarts" ) opts.restarts = std::stoul( val );
    else if ( key == "--explore-timeout" ) opts.explore_timeout = std::stoul( val );
    else if ( key == "--steps" ) opts.steps = std::stoul( val );
    else if ( key == "--steps-no-impr" ) opts.steps_no_impr = std::stoul( val );
    else if ( key == "--compress-per-step" ) opts.compress_per_step = std::stoul( val );
    else if ( key == "--spfd-k" ) opts.spfd_k = std::stoul( val );
    else if ( key == "--spfd-samples" ) opts.spfd_samples = std::stoul( val );
    else if ( key == "--spfd-beta" ) opts.spfd_beta = std::stod( val );
    else if ( key == "--spfd-max-divs" ) opts.spfd_max_divs = std::stoul( val );
    else if ( key == "--spfd-onfail" ) opts.spfd_onfail = true;
    else if ( key == "--spfd-diagnose" ) opts.spfd_diagnose = true;
    else if ( key == "--xag-flow" ) opts.xag_flow = val;
    else if ( key == "--verbose" ) g_verbose = true;
    else
    {
      std::cerr << "[mt_flow] unknown option: " << a << "\n";
      return 1;
    }
  }

  aig_network aig;
  if ( !read_input( input, aig ) )
  {
    std::cerr << "[mt_flow] could not read " << input << "\n";
    return 1;
  }

  auto const t0 = std::chrono::steady_clock::now();
  uint32_t const gates_in = aig.num_gates();
  log( fmt::format( "read {}: {} PIs, {} POs, {} AND gates", input, aig.num_pis(),
                    aig.num_pos(), gates_in ) );

  if ( !opts.explore.empty() )
  {
    uint32_t const before = aig.num_gates();
    if ( !run_explorer( aig, opts ) )
      return 1;
    aig = cleanup_dangling( aig );
    log( fmt::format( "explore/{} (cost={}): {} -> {} gates", opts.explore, opts.explore_cost,
                      before, aig.num_gates() ) );
  }

  auto const ops = split( opts.flow, ',' );
  for ( uint32_t r = 0; r < opts.rounds; ++r )
  {
    uint32_t const before_round = aig.num_gates();
    for ( auto const& op : ops )
      run_aig_op( aig, op, opts );
    aig = cleanup_dangling( aig );
    log( fmt::format( "round {}: {} -> {} gates", r, before_round, aig.num_gates() ) );
    if ( opts.rounds > 1 && aig.num_gates() >= before_round )
      break; /* converged */
  }

  uint32_t gates_out = aig.num_gates();
  klut_network klut;
  if ( opts.xag_flow.empty() )
  {
    klut = map_to_luts( aig, opts );
  }
  else
  {
    xag_network xag = cleanup_dangling<aig_network, xag_network>( aig );
    for ( auto const& xop : split( opts.xag_flow, ',' ) )
      run_xag_op( xag, xop, opts );
    gates_out = xag.num_gates();
    klut = map_to_luts( xag, opts );
  }
  depth_view<klut_network> klut_d{ klut };

  write_blif( klut, output );

  auto const secs = std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
  /* one machine-readable line on stdout so the harness log carries the numbers */
  fmt::print( "mt_flow: aig_in={} aig_out={} luts={} lut_depth={} runtime={:.2f}\n",
              gates_in, gates_out, klut.num_gates(), klut_d.depth(), secs );

  return 0;
}
