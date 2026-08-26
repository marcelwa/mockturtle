/* acd_map: decomposition-aware LUT mapping in mockturtle, for the agentic-synthesis
 * wide-cut-regeneration study.
 *
 * Usage:
 *   acd_map <input> <output.blif> [options]
 *
 *   <input>       .aig (binary AIGER) or .blif (k-LUT netlist, decomposed back to an AIG)
 *   <output.blif> a mapped BLIF with no LUT wider than the decomposition LUT size
 *
 * Options:
 *   --cut-size=N   maximum number of leaves of an enumerated cut (default 6)
 *   --acd=N        LUT size of the Ashenhurst-Curtis decomposed structure;
 *                  0 (default) disables decomposition entirely
 *   --regen        let every mapping pass generate wide (decomposable) cuts, rather
 *                  than only the first pass.  This is the intervention under study.
 *   --cut-limit=N  size of the per-node priority cut list (default 8)
 *   --delay        delay-oriented mapping (default is area-oriented)
 *   --relax=N      required-delay relaxation in % (delay-oriented only)
 *   --verbose      per-round mapper statistics on stderr
 *
 * The three arms of the experiment are
 *   ctl   --cut-size=6                 decomposition absent, objective matched
 *   acd   --cut-size=11 --acd=6        decomposition present, mechanism off
 *   wide  --cut-size=11 --acd=6 --regen
 */

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>

#include <mockturtle/algorithms/acd_expand.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
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

  std::cerr << "[acd_map] unrecognised input extension: " << path << "\n";
  return false;
}

} /* namespace */

int main( int argc, char** argv )
{
  if ( argc < 3 )
  {
    std::cerr << "usage: acd_map <input> <output.blif> [options]\n";
    return 1;
  }

  std::string const input = argv[1];
  std::string const output = argv[2];

  uint32_t cut_size = 6u;
  uint32_t acd_lut_size = 0u;
  uint32_t cut_limit = 8u;
  uint32_t relax = 0u;
  bool regen = false;
  bool area_oriented = true;
  bool verbose = false;

  for ( int i = 3; i < argc; ++i )
  {
    std::string const a = argv[i];
    auto const eq = a.find( '=' );
    std::string const key = eq == std::string::npos ? a : a.substr( 0, eq );
    std::string const val = eq == std::string::npos ? "" : a.substr( eq + 1 );

    if ( key == "--cut-size" ) cut_size = std::stoul( val );
    else if ( key == "--acd" ) acd_lut_size = std::stoul( val );
    else if ( key == "--cut-limit" ) cut_limit = std::stoul( val );
    else if ( key == "--relax" ) relax = std::stoul( val );
    else if ( key == "--regen" ) regen = true;
    else if ( key == "--delay" ) area_oriented = false;
    else if ( key == "--verbose" ) verbose = true;
    else
    {
      std::cerr << "[acd_map] unknown option: " << a << "\n";
      return 1;
    }
  }

  aig_network aig;
  if ( !read_input( input, aig ) )
  {
    std::cerr << "[acd_map] could not read " << input << "\n";
    return 1;
  }

  auto const t0 = std::chrono::steady_clock::now();

  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = cut_size;
  ps.cut_enumeration_ps.cut_limit = cut_limit;
  ps.area_oriented_mapping = area_oriented;
  ps.relax_required = relax;
  ps.acd_lut_size = acd_lut_size;
  ps.acd_regenerate_wide = regen;
  ps.verbose = verbose;

  lut_map_stats mst;
  klut_network klut = lut_map<aig_network, true>( aig, ps, &mst );

  acd_expand_stats est;
  if ( acd_lut_size != 0 )
  {
    klut = acd_expand( klut, acd_lut_size, &est );
    if ( est.num_failed != 0 )
    {
      std::cerr << fmt::format( "[acd_map] FATAL: {} wide LUTs had no decomposition\n", est.num_failed );
      return 1;
    }
  }

  /* legality guard: the harness rejects any LUT wider than the target size */
  uint32_t const limit = acd_lut_size != 0 ? acd_lut_size : cut_size;
  uint32_t widest = 0;
  klut.foreach_gate( [&]( auto const& n ) {
    widest = std::max( widest, klut.fanin_size( n ) );
  } );
  if ( widest > limit )
  {
    std::cerr << fmt::format( "[acd_map] FATAL: widest LUT has {} inputs, limit {}\n", widest, limit );
    return 1;
  }

  depth_view<klut_network> klut_d{ klut };
  write_blif( klut, output );

  auto const secs = std::chrono::duration<double>( std::chrono::steady_clock::now() - t0 ).count();
  fmt::print( "acd_map: luts={} lut_depth={} widest={} decomposed={} dec_luts={} runtime={:.2f}\n",
              klut.num_gates(), klut_d.depth(), widest, est.num_decomposed, est.num_luts_created, secs );

  return 0;
}
