/* lutmap_probe: minimal reproducer -- read an AIG verbatim, call lut_map, nothing else.
 * Used to attribute a segfault seen on `&dch -f; &put`-produced AIGs (which contain
 * unreferenced/dangling AND nodes) to lut_mapper.hpp rather than to the ACD port.
 * Throwaway.  Not part of any measurement.
 */
#include <iostream>
#include <string>

#include <fmt/format.h>
#include <lorina/aiger.hpp>

#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>

using namespace mockturtle;

int main( int argc, char** argv )
{
  if ( argc < 2 )
    return 1;
  bool const do_cleanup = argc > 2 && std::string( argv[2] ) == "--cleanup";

  aig_network aig;
  if ( lorina::read_aiger( argv[1], aiger_reader( aig ) ) != lorina::return_code::success )
    return 1;
  fmt::print( "read: size={} gates={} pis={} pos={}\n", aig.size(), aig.num_gates(),
              aig.num_pis(), aig.num_pos() );

  if ( do_cleanup )
  {
    aig = cleanup_dangling( aig );
    fmt::print( "after cleanup_dangling: gates={}\n", aig.num_gates() );
  }

  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = 6u;
  ps.cut_enumeration_ps.cut_limit = 8u;
  ps.area_oriented_mapping = true;
  auto const klut = lut_map<aig_network, true>( aig, ps );
  fmt::print( "lut_map OK: luts={}\n", klut.num_gates() );
  return 0;
}
