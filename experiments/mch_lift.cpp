/* mch_lift: cross-representation candidate generation for the agentic-synthesis project.
 *
 * BACKLOG A5 -- the "cheap version" of Hu et al., DAC'25, "Mixed Structural Choice".
 * Their operator lifts an AIG one-to-one into an XMG / MIG / XAG container, optimises
 * inside that container (so XOR and MAJ are primitives and the reachable structures
 * differ), and retains every synthesised subcircuit as a *choice* rather than replacing
 * the original. This binary does the first half only: lift, optimise, drop back to an
 * AIG, and write it out. The choice merging is then done in ABC (fraig_store /
 * fraig_restore) by the shell script that drives this.
 *
 * NOTE (honesty): the paper does not say what optimisation runs inside the XMG/MIG after
 * the one-to-one lift. The op vocabulary below is OUR design decision, not theirs.
 *
 * Usage:
 *   mch_lift <input.aig|.blif|.v> <output.aig|.blif> [options]
 *
 *   Output format follows the extension: `.aig` writes binary AIGER (the optimised AIG,
 *   unmapped, ready to be a structural choice); `.blif` runs mockturtle's LUT mapper and
 *   writes a k-LUT netlist.
 *
 * Options:
 *   --domain=aig|mig|xmg|xag   container to lift into (default aig, i.e. no lift)
 *   --flow=<a,b,c>             comma-separated ops, applied left to right
 *   --rounds=N                 repeat the flow N times (default 1); stops early on no gain
 *   --k=N                      LUT size when writing BLIF (default 6)
 *   --cut-limit=N              cut limit for the BLIF mapping (default 8)
 *   --map=area|delay           mapping style when writing BLIF (default area)
 *   --max-pis=N                resubstitution window inputs (default 8)
 *   --max-inserts=N            resubstitution insertion limit (default 2)
 *   --max-divisors=N           resubstitution divisor limit (default 150)
 *   --seed=N                   random seed where applicable
 *   --verbose                  progress on stderr
 *
 * Op vocabulary -- the SAME names in every domain, dispatched to that domain's engine.
 * This is the point of the binary: `--flow=rw,rs,ad` means "the closest thing this
 * container has to rewrite, resubstitute, algebraic depth rewriting", so the four domains
 * are run with matched effort and differ only in what the container can express.
 *
 *   rw    rewrite against the domain's 4-input NPN database
 *         aig: xag_npn (aig_complete) | mig: mig_npn | xmg: xmg3_npn | xag: xag_npn (xag_complete)
 *   rs    resubstitution (aig_resubstitution / mig_resubstitution / xmg_resubstitution /
 *         xag_resubstitution)
 *   rs2   the alternative resub engine where the domain has one (aig, mig); else = rs
 *   ad    algebraic depth rewriting (mig_/xmg_/xag_algebraic_depth_rewriting; aig: aig_balance)
 *   bal   balancing (aig_balance / xag_balance / SOP rebalancing elsewhere)
 *   sopb  SOP rebalancing at cut size k
 *   dc    don't-care based optimisation where available (xmg only); elsewhere a no-op
 *   fr    functional_reduction (aig/xag only)
 *   map   exact-library area mapping onto the domain's own gates (mig/xmg)
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
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/mig_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg3_npn.hpp>
#include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/xag_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/xag_balancing.hpp>
#include <mockturtle/algorithms/xag_resub.hpp>
#include <mockturtle/algorithms/xmg_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/xmg_optimization.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
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
  std::string domain = "aig";
  std::string flow = "rw,rs,rw";
  uint32_t rounds = 1u;
  uint32_t k = 6u;
  uint32_t cut_limit = 8u;
  std::string map_style = "area";
  uint32_t max_pis = 8u;
  uint32_t max_inserts = 2u;
  uint32_t max_divisors = 150u;
  uint32_t seed = 1u;
};

std::vector<std::string> split( std::string const& s, char sep )
{
  std::vector<std::string> out;
  std::string cur;
  std::istringstream is( s );
  while ( std::getline( is, cur, sep ) )
  {
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
    std::cerr << "[mch_lift] " << msg << "\n";
}

bool read_input( std::string const& path, aig_network& aig )
{
  if ( ends_with( path, ".aig" ) || ends_with( path, ".aag" ) )
    return lorina::read_aiger( path, aiger_reader( aig ) ) == lorina::return_code::success;
  if ( ends_with( path, ".v" ) || ends_with( path, ".verilog" ) )
    return lorina::read_verilog( path, verilog_reader( aig ) ) == lorina::return_code::success;
  if ( ends_with( path, ".blif" ) )
  {
    klut_network klut;
    names_view<klut_network> named{ klut };
    if ( lorina::read_blif( path, blif_reader( named ) ) != lorina::return_code::success )
      return false;
    aig = convert_klut_to_graph<aig_network>( named );
    return true;
  }
  std::cerr << "[mch_lift] unrecognised input extension: " << path << "\n";
  return false;
}

resubstitution_params resub_ps( options const& opts )
{
  resubstitution_params ps;
  ps.max_pis = opts.max_pis;
  ps.max_inserts = opts.max_inserts;
  ps.max_divisors = opts.max_divisors;
  return ps;
}

template<class Ntk>
void sop_rebalance( Ntk& ntk, options const& opts )
{
  sop_rebalancing<Ntk> balance_fn;
  balancing_params bps;
  bps.cut_enumeration_ps.cut_size = opts.k;
  ntk = balancing( ntk, { balance_fn }, bps );
}

/* ------------------------------------------------------------------ AIG ops */

void run_op( aig_network& aig, std::string const& op, options const& opts )
{
  if ( op == "rw" )
  {
    xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
    exact_library_params eps;
    exact_library<aig_network> lib( resyn, eps );
    rewrite_params ps;
    rewrite( aig, lib, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "rs" || op == "rs2" )
  {
    auto ps = resub_ps( opts );
    depth_view d{ aig };
    fanout_view f{ d };
    if ( op == "rs" )
      aig_resubstitution( f, ps );
    else
      aig_resubstitution2( f, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "ad" || op == "bal" )
  {
    aig_balancing_params ps;
    ps.minimize_levels = true;
    aig_balance( aig, ps );
  }
  else if ( op == "sopb" )
  {
    sop_rebalance( aig, opts );
  }
  else if ( op == "rf" )
  {
    sop_factoring<aig_network> resyn;
    refactoring_params ps;
    ps.max_pis = 10u;
    refactoring( aig, resyn, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "fr" )
  {
    functional_reduction_params ps;
    functional_reduction( aig, ps );
    aig = cleanup_dangling( aig );
  }
  else if ( op == "dc" || op == "map" )
  {
    /* no AIG equivalent; deliberately a no-op so a flow string can be shared */
  }
  else
  {
    std::cerr << "[mch_lift] unknown aig op: " << op << "\n";
  }
}

/* ------------------------------------------------------------------ MIG ops */

void run_op( mig_network& mig, std::string const& op, options const& opts )
{
  if ( op == "rw" || op == "map" )
  {
    mig_npn_resynthesis resyn{ true };
    exact_library_params eps;
    exact_library<mig_network> lib( resyn, eps );
    if ( op == "rw" )
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
  else if ( op == "rs" || op == "rs2" )
  {
    auto ps = resub_ps( opts );
    depth_view d{ mig };
    fanout_view f{ d };
    if ( op == "rs" )
      mig_resubstitution( f, ps );
    else
      mig_resubstitution2( f, ps );
    mig = cleanup_dangling( mig );
  }
  else if ( op == "ad" )
  {
    depth_view d{ mig };
    mig_algebraic_depth_rewriting( d );
    mig = cleanup_dangling( mig );
  }
  else if ( op == "bal" || op == "sopb" )
  {
    sop_rebalance( mig, opts );
  }
  else if ( op == "dc" || op == "fr" )
  {
    /* no MIG equivalent */
  }
  else
  {
    std::cerr << "[mch_lift] unknown mig op: " << op << "\n";
  }
}

/* ------------------------------------------------------------------ XMG ops */

void run_op( xmg_network& xmg, std::string const& op, options const& opts )
{
  if ( op == "rw" || op == "map" )
  {
    xmg3_npn_resynthesis<xmg_network> resyn;
    exact_library_params eps;
    eps.np_classification = false;
    exact_library<xmg_network> lib( resyn, eps );
    if ( op == "rw" )
    {
      rewrite_params ps;
      rewrite( xmg, lib, ps );
      xmg = cleanup_dangling( xmg );
    }
    else
    {
      map_params mps;
      mps.skip_delay_round = true;
      mps.required_time = std::numeric_limits<double>::max();
      xmg = map( xmg, lib, mps );
    }
  }
  else if ( op == "rs" || op == "rs2" )
  {
    auto ps = resub_ps( opts );
    depth_view d{ xmg };
    fanout_view f{ d };
    xmg_resubstitution( f, ps );
    xmg = cleanup_dangling( xmg );
  }
  else if ( op == "ad" )
  {
    depth_view d{ xmg };
    xmg_algebraic_depth_rewriting( d );
    xmg = cleanup_dangling( xmg );
  }
  else if ( op == "dc" )
  {
    xmg = xmg_dont_cares_optimization( xmg );
  }
  else if ( op == "bal" || op == "sopb" )
  {
    sop_rebalance( xmg, opts );
  }
  else if ( op == "fr" )
  {
    /* no XMG equivalent */
  }
  else
  {
    std::cerr << "[mch_lift] unknown xmg op: " << op << "\n";
  }
}

/* ------------------------------------------------------------------ XAG ops */

void run_op( xag_network& xag, std::string const& op, options const& opts )
{
  if ( op == "rw" || op == "map" )
  {
    xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete> resyn;
    exact_library_params eps;
    exact_library<xag_network> lib( resyn, eps );
    if ( op == "rw" )
    {
      rewrite_params ps;
      rewrite( xag, lib, ps );
      xag = cleanup_dangling( xag );
    }
    else
    {
      map_params mps;
      mps.skip_delay_round = true;
      mps.required_time = std::numeric_limits<double>::max();
      xag = map( xag, lib, mps );
    }
  }
  else if ( op == "rs" || op == "rs2" )
  {
    auto ps = resub_ps( opts );
    depth_view d{ xag };
    fanout_view f{ d };
    xag_resubstitution( f, ps );
    xag = cleanup_dangling( xag );
  }
  else if ( op == "ad" )
  {
    depth_view d{ xag };
    xag_algebraic_depth_rewriting( d );
    xag = cleanup_dangling( xag );
  }
  else if ( op == "bal" )
  {
    xag_balancing_params ps;
    xag_balance( xag, ps );
    xag = cleanup_dangling( xag );
  }
  else if ( op == "sopb" )
  {
    sop_rebalance( xag, opts );
  }
  else if ( op == "fr" )
  {
    functional_reduction_params ps;
    functional_reduction( xag, ps );
    xag = cleanup_dangling( xag );
  }
  else if ( op == "dc" )
  {
    /* no XAG equivalent */
  }
  else
  {
    std::cerr << "[mch_lift] unknown xag op: " << op << "\n";
  }
}

/* --------------------------------------------------------------- excursion */

/* Lift the AIG one-to-one into `Ntk`, run the flow there, drop back to an AIG.
 * `cleanup_dangling<A,B>` is exactly the one-to-one container change the paper
 * describes: it walks the source in topological order and re-creates each AND as an
 * AND in the destination. No resynthesis happens in the lift itself. */
template<class Ntk>
void excursion( aig_network& aig, options const& opts )
{
  Ntk ntk = cleanup_dangling<aig_network, Ntk>( aig );
  log( fmt::format( "lifted to {} gates", ntk.num_gates() ) );

  auto const ops = split( opts.flow, ',' );
  for ( uint32_t r = 0; r < opts.rounds; ++r )
  {
    uint32_t const before_round = ntk.num_gates();
    for ( auto const& op : ops )
    {
      auto const before = ntk.num_gates();
      run_op( ntk, op, opts );
      log( fmt::format( "  {}: {} -> {}", op, before, ntk.num_gates() ) );
    }
    ntk = cleanup_dangling( ntk );
    if ( opts.rounds > 1 && ntk.num_gates() >= before_round )
      break;
  }

  aig = cleanup_dangling<Ntk, aig_network>( ntk );
}

klut_network map_to_luts( aig_network const& aig, options const& opts )
{
  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = opts.k;
  ps.cut_enumeration_ps.cut_limit = opts.cut_limit;
  ps.recompute_cuts = true;
  ps.cut_expansion = true;
  ps.area_oriented_mapping = ( opts.map_style == "area" );
  return lut_map( aig, ps );
}

} // namespace

int main( int argc, char** argv )
{
  if ( argc < 3 )
  {
    std::cerr << "usage: mch_lift <input> <output.aig|.blif> [--domain=aig|mig|xmg|xag] "
                 "[--flow=...] [--rounds=N] [--k=N] [--cut-limit=N] [--map=area|delay] "
                 "[--max-pis=N] [--max-inserts=N] [--max-divisors=N] [--seed=N] [--verbose]\n";
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

    if ( key == "--domain" ) opts.domain = val;
    else if ( key == "--flow" ) opts.flow = val;
    else if ( key == "--rounds" ) opts.rounds = std::stoul( val );
    else if ( key == "--k" ) opts.k = std::stoul( val );
    else if ( key == "--cut-limit" ) opts.cut_limit = std::stoul( val );
    else if ( key == "--map" ) opts.map_style = val;
    else if ( key == "--max-pis" ) opts.max_pis = std::stoul( val );
    else if ( key == "--max-inserts" ) opts.max_inserts = std::stoul( val );
    else if ( key == "--max-divisors" ) opts.max_divisors = std::stoul( val );
    else if ( key == "--seed" ) opts.seed = std::stoul( val );
    else if ( key == "--verbose" ) g_verbose = true;
    else
    {
      std::cerr << "[mch_lift] unknown option: " << a << "\n";
      return 1;
    }
  }

  aig_network aig;
  if ( !read_input( input, aig ) )
  {
    std::cerr << "[mch_lift] could not read " << input << "\n";
    return 1;
  }

  auto const t0 = std::chrono::steady_clock::now();
  uint32_t const gates_in = aig.num_gates();
  log( fmt::format( "read {}: {} PIs, {} POs, {} AND gates", input, aig.num_pis(),
                    aig.num_pos(), gates_in ) );

  if ( opts.domain == "aig" )
    excursion<aig_network>( aig, opts );
  else if ( opts.domain == "mig" )
    excursion<mig_network>( aig, opts );
  else if ( opts.domain == "xmg" )
    excursion<xmg_network>( aig, opts );
  else if ( opts.domain == "xag" )
    excursion<xag_network>( aig, opts );
  else
  {
    std::cerr << "[mch_lift] unknown --domain: " << opts.domain << "\n";
    return 1;
  }

  aig = cleanup_dangling( aig );

  uint32_t luts = 0u, lut_depth = 0u;
  if ( ends_with( output, ".blif" ) )
  {
    auto const klut = map_to_luts( aig, opts );
    depth_view<klut_network> klut_d{ klut };
    luts = klut.num_gates();
    lut_depth = klut_d.depth();
    write_blif( klut, output );
  }
  else
  {
    write_aiger( aig, output );
  }

  depth_view<aig_network> aig_d{ aig };
  auto const secs = std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
  fmt::print( "mch_lift: domain={} aig_in={} aig_out={} aig_depth={} luts={} lut_depth={} runtime={:.2f}\n",
              opts.domain, gates_in, aig.num_gates(), aig_d.depth(), luts, lut_depth, secs );

  return 0;
}
