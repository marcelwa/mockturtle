/* acd_resyn2: acd_resyn.cpp plus ONE additional engine, `dsdacd`.
 *
 * Kept as a separate file, and therefore a separate binary, on purpose: the confirmatory
 * evaluation-A/B sweeps of 2026-08-26-p2-acd-resyn were run against `acd_resyn` and that
 * binary must stay byte-identical to what its ledger rows record.  This one exists only for
 * the POST-HOC, exploratory follow-up arm, and nothing scored with it may be reported as
 * part of the pre-registered result.
 *
 * `dsdacd` is `dsd_resynthesis` with `acd_resynthesis` as its prime-remainder handler,
 * which in turn falls back to the same Shannon+NPN engine as everything else.  It is
 * motivated by an observation from the confirmatory run and by nothing else: ACD beat DSD
 * on the median benchmark but lost by 1.89x on `div` and 1.63x on `sqrt`, which are exactly
 * the two benchmarks where DSD is at its best against the shared fallback (0.39x and
 * 0.59x).  The reading is that ACD is a generalisation in EXPRESSIVENESS but not in OUTPUT:
 * having committed to a free-set/bound-set split it does not recover the full disjoint
 * decomposition when one exists.  Composing them tests exactly that reading.
 *
 * Original header follows.
 *
 * acd_resyn: a choice of resynthesis engine applied to a logic network, so that several
 * engines can be A/B-ed as drop-in replacements for one another on identical input.
 *
 * The head-to-head this exists for is `acd_resynthesis` (node_resynthesis/acd.hpp) against
 * the engine it generalises, `dsd_resynthesis`.  Every engine below is given the *same*
 * fallback (`shannon_resynthesis` down to 4 variables, then a COMPLETE 4-input NPN
 * database), so the difference between two arms is the decomposition and nothing else.
 * The same object also realises the blocks of the ACD cascade.
 *
 * Usage:
 *   acd_resyn <input> <output> [options]
 *
 *   <input>    .aig (binary AIGER) or .blif (k-LUT netlist, decomposed back to an AIG)
 *   <output>   .aig (default) or .blif, per --emit
 *
 * Options:
 *   --engine=E     direct | shannon | dsd | acd | dsdacd | bidec | sopf   (default acd)
 *   --form=F       noderesyn (default) | cutrw
 *   --map-k=K      LUT size of the intermediate k-LUT network (noderesyn), or the cut size
 *                  (cutrw).  Default 11 -- the ACD operating point; ACD accepts nothing
 *                  above 11 variables (2026-08-26-acd-resyn-census).
 *   --cut-limit=N  priority-cut list size (cutrw only, default 12)
 *   --allow-zero-gain   accept zero-gain candidates (cutrw only)
 *   --acd-lut-size=N    block size of the ACD cascade (default 6, the census's
 *                       operating point; the blocks are realised by the shared fallback)
 *   --emit=aig|blif     what to write (default aig)
 *   --emit-k=K          LUT size when --emit=blif (default 6)
 *   --verbose
 *
 * `--form=noderesyn` is the primary form:
 *
 *     AIG -> lut_map(K=--map-k, area) -> k-LUT network -> node_resynthesis<aig> -> AIG
 *
 * The mapping step is what "asks for wide cuts deliberately"; mockturtle's stock
 * smallest-first cut order offers a >6-variable cut at only 17.2 % of nodes.
 *
 * `--form=cutrw` runs `cut_rewriting` directly on the input AIG instead.  Note that
 * `cut_rewriting`'s acceptance gain is measured in **AIG nodes** while ACD trades AIG nodes
 * for 6-LUTs, so this form is expected to be near-null; it is measured, not fixed.
 */

#include <chrono>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>

/* MUST come first: the vendored ACD sources are namespaced by an ABC macro that
 * `acd_namespace.hpp` defines only if nothing else has claimed it yet.  ABC's own
 * `abc_namespaces.h` -- pulled in transitively by cut_rewriting/sop_factoring via the SAT and
 * ESOP libraries -- redefines it to `pabc`, so the whole ACD chain must be included before
 * any of them or the vendored decomposer lands in `pabc::acd` and stops compiling. */
#include <mockturtle/algorithms/node_resynthesis/acd.hpp>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/bidecomposition.hpp>
#include <mockturtle/algorithms/node_resynthesis/dsd.hpp>
#include <mockturtle/algorithms/node_resynthesis/shannon.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>

using namespace mockturtle;

namespace
{

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

  std::cerr << "[acd_resyn] unrecognised input extension: " << path << "\n";
  return false;
}

/* the shared fallback: Shannon down to 4 variables, then a COMPLETE 4-input NPN database.
 * The default database is `xag_incomplete` and misses roughly two thirds of the classes,
 * failing silently -- that already cost the Phase 1 correctness gate a day. */
using npn_t = xag_npn_resynthesis<aig_network, xag_network, xag_npn_db_kind::xag_complete>;
using fallback_t = shannon_resynthesis<aig_network, npn_t>;

double secs_since( std::chrono::steady_clock::time_point t )
{
  return std::chrono::duration<double>( std::chrono::steady_clock::now() - t ).count();
}

} /* namespace */

int main( int argc, char** argv )
{
  if ( argc < 3 )
  {
    std::cerr << "usage: acd_resyn2 <input> <output> [options]\n";
    return 1;
  }

  std::string const input = argv[1];
  std::string const output = argv[2];

  std::string engine = "acd";
  std::string form = "noderesyn";
  std::string emit = "aig";
  uint32_t map_k = 11u;
  uint32_t cut_limit = 12u;
  uint32_t acd_lut_size = 6u;
  uint32_t emit_k = 6u;
  bool allow_zero_gain = false;
  bool verbose = false;

  for ( int i = 3; i < argc; ++i )
  {
    std::string const a = argv[i];
    auto const eq = a.find( '=' );
    std::string const key = eq == std::string::npos ? a : a.substr( 0, eq );
    std::string const val = eq == std::string::npos ? "" : a.substr( eq + 1 );

    if ( key == "--engine" ) engine = val;
    else if ( key == "--form" ) form = val;
    else if ( key == "--emit" ) emit = val;
    else if ( key == "--map-k" ) map_k = std::stoul( val );
    else if ( key == "--cut-limit" ) cut_limit = std::stoul( val );
    else if ( key == "--acd-lut-size" ) acd_lut_size = std::stoul( val );
    else if ( key == "--emit-k" ) emit_k = std::stoul( val );
    else if ( key == "--allow-zero-gain" ) allow_zero_gain = true;
    else if ( key == "--verbose" ) verbose = true;
    else
    {
      std::cerr << "[acd_resyn] unknown option: " << a << "\n";
      return 1;
    }
  }

  if ( engine != "direct" && engine != "shannon" && engine != "dsd" && engine != "acd" &&
       engine != "bidec" && engine != "sopf" && engine != "dsdacd" )
  {
    std::cerr << "[acd_resyn] unknown engine: " << engine << "\n";
    return 1;
  }
  if ( form != "noderesyn" && form != "cutrw" )
  {
    std::cerr << "[acd_resyn] unknown form: " << form << "\n";
    return 1;
  }
  if ( acd_lut_size < 2u || acd_lut_size > 6u )
  {
    std::cerr << "[acd_resyn] --acd-lut-size must be in [2,6]\n";
    return 1;
  }
  if ( emit != "aig" && emit != "blif" )
  {
    std::cerr << "[acd_resyn] unknown emit mode: " << emit << "\n";
    return 1;
  }

  aig_network aig_in;
  if ( !read_input( input, aig_in ) )
  {
    std::cerr << "[acd_resyn] could not read " << input << "\n";
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

  /* the engines.  All of them share `fb`, so an arm differs from another arm only in the
   * decomposition it puts in front of the shared fallback. */
  npn_t npn;
  /* The Shannon threshold is the NPN DATABASE's width (4), never the ACD block size.
   * `xag_npn_resynthesis` reads past its 4-input database rather than asserting, so handing
   * it a 5- or 6-variable block is a segfault, not an error.  With the threshold at 4 the
   * fallback is total for any arity and every leaf call is exactly 4 variables wide. */
  fallback_t fb{ 4u, &npn };

  acd_resynthesis_params acd_ps;
  acd_ps.lut_size = acd_lut_size;
  acd_ps.max_num_vars = 11u;
  acd_ps.use_fallback = true; /* node_resynthesis needs total coverage */
  acd_resynthesis_stats acd_st;

  /* `cut_rewriting` wants a *candidate-offering* operator, so no fallback there: a decline
   * must emit nothing and leave the original structure in place. */
  acd_resynthesis_params acd_ps_cr = acd_ps;
  acd_ps_cr.use_fallback = false;

  double secs_map = 0.0;
  double secs_resyn = 0.0;
  uint32_t klut_size = 0u;
  uint32_t klut_max_fanin = 0u;

  aig_network aig_out;

  try
  {
    if ( form == "noderesyn" )
    {
      auto const t_map = std::chrono::steady_clock::now();
      lut_map_params mps;
      mps.cut_enumeration_ps.cut_size = map_k;
      mps.cut_enumeration_ps.cut_limit = cut_limit;
      mps.area_oriented_mapping = true;
      mps.verbose = verbose;
      klut_network klut = lut_map<aig_network, true>( aig_in, mps );
      secs_map = secs_since( t_map );

      klut_size = klut.num_gates();
      klut.foreach_gate( [&]( auto const& n ) {
        klut_max_fanin = std::max( klut_max_fanin, klut.fanin_size( n ) );
      } );

      auto const t_resyn = std::chrono::steady_clock::now();
      if ( engine == "direct" )
      {
        aig_out = convert_klut_to_graph<aig_network>( klut );
      }
      else if ( engine == "shannon" )
      {
        aig_out = node_resynthesis<aig_network>( klut, fb );
      }
      else if ( engine == "dsd" )
      {
        dsd_resynthesis<aig_network, fallback_t> resyn( fb );
        aig_out = node_resynthesis<aig_network>( klut, resyn );
      }
      else if ( engine == "acd" )
      {
        acd_resynthesis<aig_network, fallback_t> resyn( fb, acd_ps, &acd_st );
        aig_out = node_resynthesis<aig_network>( klut, resyn );
      }
      else if ( engine == "dsdacd" )
      {
        /* DSD first; whatever it cannot peel disjointly is handed to ACD, and only what ACD
         * also declines reaches the shared Shannon+NPN fallback. */
        acd_resynthesis<aig_network, fallback_t> inner( fb, acd_ps, &acd_st );
        dsd_resynthesis<aig_network, decltype( inner )> resyn( inner );
        aig_out = node_resynthesis<aig_network>( klut, resyn );
      }
      else if ( engine == "bidec" )
      {
        bidecomposition_resynthesis<aig_network> resyn;
        aig_out = node_resynthesis<aig_network>( klut, resyn );
      }
      else /* sopf */
      {
        sop_factoring<aig_network> resyn;
        aig_out = node_resynthesis<aig_network>( klut, resyn );
      }
      secs_resyn = secs_since( t_resyn );
    }
    else /* cutrw */
    {
      cut_rewriting_params cps;
      cps.cut_enumeration_ps.cut_size = map_k;
      cps.cut_enumeration_ps.cut_limit = cut_limit;
      cps.cut_enumeration_ps.minimize_truth_table = true;
      cps.allow_zero_gain = allow_zero_gain;
      cps.verbose = verbose;

      auto const t_resyn = std::chrono::steady_clock::now();
      if ( engine == "direct" )
      {
        aig_out = aig_in;
      }
      else if ( engine == "shannon" )
      {
        aig_out = cut_rewriting( aig_in, fb, cps );
      }
      else if ( engine == "dsd" )
      {
        dsd_resynthesis<aig_network, fallback_t> resyn( fb );
        aig_out = cut_rewriting( aig_in, resyn, cps );
      }
      else if ( engine == "acd" )
      {
        acd_resynthesis<aig_network, fallback_t> resyn( fb, acd_ps_cr, &acd_st );
        aig_out = cut_rewriting( aig_in, resyn, cps );
      }
      else if ( engine == "dsdacd" )
      {
        acd_resynthesis<aig_network, fallback_t> inner( fb, acd_ps, &acd_st );
        dsd_resynthesis<aig_network, decltype( inner )> resyn( inner );
        aig_out = cut_rewriting( aig_in, resyn, cps );
      }
      else if ( engine == "bidec" )
      {
        bidecomposition_resynthesis<aig_network> resyn;
        aig_out = cut_rewriting( aig_in, resyn, cps );
      }
      else /* sopf */
      {
        sop_factoring<aig_network> resyn;
        aig_out = cut_rewriting( aig_in, resyn, cps );
      }
      secs_resyn = secs_since( t_resyn );
    }
  }
  catch ( std::exception const& e )
  {
    std::cerr << "[acd_resyn] FATAL: engine " << engine << " threw: " << e.what() << "\n";
    return 1;
  }

  aig_out = cleanup_dangling( aig_out );

  if ( acd_st.num_block_failures != 0u )
  {
    std::cerr << fmt::format( "[acd_resyn] FATAL: {} blocks could not be realised (first at {} fanins)\n",
                              acd_st.num_block_failures, acd_st.first_failed_fanins );
    return 1;
  }

  if ( aig_out.num_pis() != pi_in || aig_out.num_pos() != po_in )
  {
    std::cerr << fmt::format( "[acd_resyn] FATAL: interface changed, PI {}->{} PO {}->{}\n",
                              pi_in, aig_out.num_pis(), po_in, aig_out.num_pos() );
    return 1;
  }

  uint32_t aig_out_size = aig_out.num_gates();
  uint32_t aig_out_depth = 0;
  {
    depth_view<aig_network> d{ aig_out };
    aig_out_depth = d.depth();
  }

  if ( emit == "aig" )
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
      std::cerr << fmt::format( "[acd_resyn] FATAL: widest LUT has {} inputs, limit {}\n", widest, emit_k );
      return 1;
    }
    write_blif( mapped, output );
  }

  double const secs_total = secs_since( t_start );

  fmt::print( stderr,
              "[acd_resyn] engine={} form={} map_k={} pi={} po={} "
              "aig_in_size={} aig_in_depth={} klut_size={} klut_max_fanin={} "
              "aig_out_size={} aig_out_depth={} "
              "acd_accepted={} acd_declined={} acd_trivial={} acd_blocks={} acd_block_failures={} "
              "acd_max_block_fanins={} "
              "secs_map={:.3f} secs_resyn={:.3f} secs_total={:.3f}\n",
              engine, form, map_k, pi_in, po_in,
              aig_in_size, aig_in_depth, klut_size, klut_max_fanin,
              aig_out_size, aig_out_depth,
              acd_st.num_accepted, acd_st.num_declined, acd_st.num_trivial, acd_st.num_blocks,
              acd_st.num_block_failures, acd_st.max_block_fanins,
              secs_map, secs_resyn, secs_total );

  return 0;
}
