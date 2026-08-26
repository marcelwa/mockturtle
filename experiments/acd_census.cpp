/* acd_census: how big is the surface for ACD as an AIG resynthesis operator?
 *
 * Exploratory instrumentation, NOT a scored experiment.  It produces no netlist,
 * runs no CEC, and nothing it prints is a QoR result.  Its single job is to
 * answer, before an engine is written:
 *
 *     on a real AIG, what fraction of nodes have a k-cut whose function the
 *     Ashenhurst-Curtis decomposer accepts, as k runs from kmin to kmax?
 *
 * This is the analogue, for the resynthesis setting, of the mapper-side cover
 * census of 2026-08-26 (experiments/2026-08-26-f1-acd-encoding-headroom), and
 * it exists because the standing instruction in the project record is "any
 * future ACD arm should run the cover census first".
 *
 * Usage:
 *   acd_census <input.aig|.blif> [options]
 *
 * Options:
 *   --name=S        label for the CSV rows (default: the input path)
 *   --kmin=N        smallest cut size to enumerate (default 7)
 *   --kmax=N        largest  cut size to enumerate (default 16, the compile-time
 *                   maximum of mockturtle's cut_enumeration)
 *   --cut-limit=N   per-node priority cut list size (default 8)
 *   --lut-size=N    target LUT size handed to ACD (default 6)
 *   --csv=PATH      append machine-readable rows here
 *   --hist          print the ACD cost / support histograms per k
 *   --wide-first    keep the WIDEST cuts in each node's priority list instead of
 *                   the smallest (mockturtle's default).  Without this the
 *                   enumerator never produces cuts above ~10 leaves and the
 *                   k >= 11 part of the curve is an artefact of the cut order.
 *
 * What is measured, per k
 *   offered      cuts enumerated at this k over all AND nodes, excluding the
 *                trivial self-cut
 *   wide         of those, the ones whose *support* exceeds lut_size, i.e. the
 *                ones for which a decomposition is needed at all.  ACD is only
 *                asked about these; a cut with support <= lut_size is one LUT
 *                already and is not part of the surface.
 *   accepted     wide cuts for which acd_iface::evaluate returned >= 0
 *   improving    accepted cuts whose LUT count is strictly below the number of
 *                6-LUT cells the baseline cover currently uses inside the same
 *                window (root node down to the cut leaves)
 *
 *   node_any     nodes with at least one accepted cut
 *   node_imp     nodes with at least one improving accepted cut
 *
 * The baseline cover comes from mockturtle's own lut_mapping at K=6, area
 * oriented, and "cells inside the window" counts cover cell roots strictly
 * between the cut leaves and the node, plus the node itself.  That is a proxy,
 * not an exact accounting -- a cell rooted inside the window may reach below the
 * cut leaves -- and it is quoted only as an order-of-magnitude comparator.
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include <mockturtle/algorithms/acd/acd_wrapper.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/views/names_view.hpp>

namespace mockturtle
{

/*! \brief Cut data that makes the priority list keep the WIDEST cuts.
 *
 * mockturtle's default cut order is smallest-first (utils/cuts.hpp), which is
 * right for a mapper and wrong for this census: with a 25-slot priority list
 * every slot goes to a 2-6 leaf cut, no wide cut survives to be merged upwards,
 * and the enumerator stops producing cuts above ~10 leaves entirely.  That would
 * make the k >= 11 part of the acceptance curve an artefact of the enumerator
 * rather than a statement about ACD.  Reversing the order measures the other
 * extreme; both are reported.
 */
struct census_wide_cut
{
};

template<bool ComputeTruth>
bool operator<( cut_type<ComputeTruth, census_wide_cut> const& c1,
                cut_type<ComputeTruth, census_wide_cut> const& c2 )
{
  return c1.size() > c2.size();
}

} /* namespace mockturtle */

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

  std::cerr << "[acd_census] unrecognised input extension: " << path << "\n";
  return false;
}

/* Per-k accumulators. */
struct k_stats
{
  uint64_t offered{ 0 };   /* non-trivial cuts enumerated                        */
  uint64_t wide{ 0 };      /* support > lut_size, i.e. ACD was actually asked     */
  uint64_t over_ceiling{ 0 }; /* wide cuts whose support > ACD's 11-var ceiling   */
  uint64_t accepted{ 0 };
  uint64_t improving{ 0 };
  uint32_t node_wide{ 0 }; /* nodes offered at least one wide cut                 */
  uint32_t node_any{ 0 };  /* nodes with at least one accepted cut                */
  uint32_t node_imp{ 0 };  /* nodes with at least one improving accepted cut      */
  uint64_t cost_sum{ 0 };
  uint64_t base_sum{ 0 };  /* baseline cells for the same accepted windows        */
  uint64_t cone_sum{ 0 };  /* AIG nodes inside the same accepted windows          */
  std::map<uint32_t, uint64_t> cost_hist;      /* ACD LUT count -> count          */
  std::map<uint32_t, uint64_t> support_wide;   /* support size -> wide cuts       */
  std::map<uint32_t, uint64_t> support_acc;    /* support size -> accepted cuts   */
  std::map<uint32_t, uint64_t> leaves_hist;    /* cut leaf count -> cuts offered  */
};

} /* namespace */

int main( int argc, char** argv )
{
  if ( argc < 2 )
  {
    std::cerr << "usage: acd_census <input.aig|.blif> [options]\n";
    return 1;
  }

  std::string const input = argv[1];
  std::string name = input;
  uint32_t kmin = 7u, kmax = 16u, cut_limit = 8u, lut_size = 6u;
  std::string csv;
  bool hist = false;
  bool wide_first = false;

  for ( int i = 2; i < argc; ++i )
  {
    std::string const a = argv[i];
    auto const eq = a.find( '=' );
    std::string const key = eq == std::string::npos ? a : a.substr( 0, eq );
    std::string const val = eq == std::string::npos ? "" : a.substr( eq + 1 );

    if ( key == "--name" ) name = val;
    else if ( key == "--kmin" ) kmin = std::stoul( val );
    else if ( key == "--kmax" ) kmax = std::stoul( val );
    else if ( key == "--cut-limit" ) cut_limit = std::stoul( val );
    else if ( key == "--lut-size" ) lut_size = std::stoul( val );
    else if ( key == "--csv" ) csv = val;
    else if ( key == "--hist" ) hist = true;
    else if ( key == "--wide-first" ) wide_first = true;
    else
    {
      std::cerr << "[acd_census] unknown option: " << a << "\n";
      return 1;
    }
  }

  if ( kmax > 16u )
  {
    std::cerr << "[acd_census] kmax capped at 16 (cut_enumeration's max_cut_size)\n";
    kmax = 16u;
  }

  aig_network aig;
  if ( !read_input( input, aig ) )
  {
    std::cerr << "[acd_census] could not read " << input << "\n";
    return 1;
  }

  uint32_t const num_gates = aig.num_gates();

  /* -------- baseline 6-LUT cover, used only as the local-cost comparator ------ */
  std::vector<uint8_t> is_cell( aig.size(), 0u );
  uint32_t baseline_cells = 0;
  {
    mapping_view<aig_network, false> mapped{ aig };
    lut_mapping_params mps;
    mps.cut_enumeration_ps.cut_size = lut_size;
    mps.cut_enumeration_ps.cut_limit = 8u;
    lut_mapping<mapping_view<aig_network, false>, false>( mapped, mps );
    aig.foreach_gate( [&]( auto const& n ) {
      if ( mapped.is_cell_root( n ) )
      {
        is_cell[aig.node_to_index( n )] = 1u;
        ++baseline_cells;
      }
    } );
  }

  fmt::print( "# acd_census {}  pis={} pos={} and={} baseline6lut={} order={} cut_limit={}\n",
              name, aig.num_pis(), aig.num_pos(), num_gates, baseline_cells,
              wide_first ? "wide-first" : "small-first", cut_limit );
  fmt::print( "# NOT A RESULT: exploratory instrumentation, no CEC, no netlist.\n" );

  /* re-usable window walker: counts cover cells and AIG nodes between n and a cut */
  std::vector<uint32_t> mark( aig.size(), 0u );
  uint32_t gen = 0u;
  std::vector<aig_network::node> stack;

  /* returns { cover cells rooted in the window, AIG nodes in the window } */
  auto window_cost = [&]( aig_network::node const& root, auto const& cut ) -> std::pair<uint32_t, uint32_t> {
    ++gen;
    for ( auto l : cut )
      mark[l] = gen;
    uint32_t cells = 0, nodes = 0;
    stack.clear();
    stack.push_back( root );
    mark[aig.node_to_index( root )] = gen; /* guard against re-push; counted below */
    ++nodes;
    if ( is_cell[aig.node_to_index( root )] )
      ++cells;
    for ( size_t s = 0; s < stack.size(); ++s )
    {
      auto const n = stack[s];
      aig.foreach_fanin( n, [&]( auto const& f ) {
        auto const c = aig.get_node( f );
        auto const ci = aig.node_to_index( c );
        if ( mark[ci] == gen || aig.is_constant( c ) || aig.is_ci( c ) )
          return;
        mark[ci] = gen;
        ++nodes;
        if ( is_cell[ci] )
          ++cells;
        stack.push_back( c );
      } );
    }
    return { cells, nodes };
  };

  std::ofstream csv_out;
  if ( !csv.empty() )
  {
    bool const fresh = !std::ifstream( csv ).good();
    csv_out.open( csv, std::ios::app );
    if ( fresh )
      csv_out << "bench,order,k,cut_limit,nodes,offered,wide,over_ceiling,accepted,improving,"
                 "node_wide,node_any,node_imp,acd_cost_sum,base_cell_sum,cone_aig_sum,secs\n";
  }

  fmt::print( "{:>3} {:>10} {:>10} {:>8} {:>7} | {:>7} {:>7} {:>7} | {:>7} {:>7} {:>7}\n",
              "k", "offered", "wide", "accept", "acc/wide", "nd_wide", "nd_any", "nd_imp",
              "avgcost", "avgbase", "avgcone" );

  for ( uint32_t k = kmin; k <= kmax; ++k )
  {
    auto const t0 = std::chrono::steady_clock::now();
    k_stats ks;

    cut_enumeration_params cps;
    cps.cut_size = k;
    cps.cut_limit = cut_limit;
    cps.minimize_truth_table = false;

    auto body = [&]( auto const& cuts ) {
    aig.foreach_gate( [&]( auto const& n ) {
      auto const idx = aig.node_to_index( n );
      bool node_wide = false, node_any = false, node_imp = false;

      for ( auto const* cut : cuts.cuts( idx ) )
      {
        if ( cut->size() <= 1u )
          continue; /* trivial self-cut */
        ++ks.offered;
        ++ks.leaves_hist[static_cast<uint32_t>( cut->size() )];

        auto tt = cuts.truth_table( *cut );
        kitty::min_base_inplace( tt );
        uint32_t support = 0;
        for ( uint32_t v = 0; v < tt.num_vars(); ++v )
          if ( kitty::has_var( tt, v ) )
            ++support;
        if ( support <= lut_size )
          continue; /* already one LUT: not part of the decomposition surface */

        ++ks.wide;
        ++ks.support_wide[support];
        node_wide = true;

        /* ACD's two-level engine hard-rejects above 11 variables; record that
         * separately so the k = 11 boundary is visible rather than inferred. */
        if ( support > 11u )
          ++ks.over_ceiling;

        auto small = kitty::shrink_to( tt, support < 6u ? 6u : support );
        uint32_t cost = 0;
        int const levels = acd_iface::evaluate( small._bits.data(), support, lut_size, &cost );
        if ( levels < 0 )
          continue;

        ++ks.accepted;
        ++ks.support_acc[support];
        ++ks.cost_hist[cost];
        ks.cost_sum += cost;
        node_any = true;

        auto const [base, cone] = window_cost( n, *cut );
        ks.base_sum += base;
        ks.cone_sum += cone;
        if ( cost < base )
        {
          ++ks.improving;
          node_imp = true;
        }
      }

      if ( node_wide ) ++ks.node_wide;
      if ( node_any )  ++ks.node_any;
      if ( node_imp )  ++ks.node_imp;
    } );
    };

    if ( wide_first )
      body( cut_enumeration<aig_network, true, census_wide_cut>( aig, cps ) );
    else
      body( cut_enumeration<aig_network, true>( aig, cps ) );

    auto const secs = std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();

    double const acc_rate = ks.wide ? 100.0 * double( ks.accepted ) / double( ks.wide ) : 0.0;
    double const avgcost = ks.accepted ? double( ks.cost_sum ) / double( ks.accepted ) : 0.0;
    double const avgbase = ks.accepted ? double( ks.base_sum ) / double( ks.accepted ) : 0.0;
    double const avgcone = ks.accepted ? double( ks.cone_sum ) / double( ks.accepted ) : 0.0;

    fmt::print( "{:>3} {:>10} {:>10} {:>8} {:>6.1f}% | {:>7} {:>7} {:>7} | {:>7.2f} {:>7.2f} {:>7.2f}  ({:.1f}s)\n",
                k, ks.offered, ks.wide, ks.accepted, acc_rate,
                ks.node_wide, ks.node_any, ks.node_imp, avgcost, avgbase, avgcone, secs );

    if ( hist )
    {
      fmt::print( "      leaves offered:" );
      for ( auto const& [s, c] : ks.leaves_hist )
        fmt::print( " {}:{}", s, c );
      fmt::print( "\n      support(wide/accepted):" );
      for ( auto const& [s, c] : ks.support_wide )
        fmt::print( " {}:{}/{}", s, c, ks.support_acc.count( s ) ? ks.support_acc.at( s ) : 0 );
      fmt::print( "\n      acd_cost hist:" );
      for ( auto const& [c, n] : ks.cost_hist )
        fmt::print( " {}:{}", c, n );
      fmt::print( "\n" );
    }

    if ( csv_out.is_open() )
      csv_out << fmt::format( "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{:.2f}\n",
                              name, wide_first ? "wide-first" : "small-first", k, cut_limit, num_gates, ks.offered, ks.wide,
                              ks.over_ceiling, ks.accepted, ks.improving,
                              ks.node_wide, ks.node_any, ks.node_imp,
                              ks.cost_sum, ks.base_sum, ks.cone_sum, secs );
  }

  return 0;
}
