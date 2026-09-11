#include <mockturtle/algorithms/exorcism.hpp>

static_assert( sizeof( abc::exorcism::ABC_PTRINT_T ) == sizeof( void* ) );

bool esop_check()
{
  kitty::dynamic_truth_table function( 4u );
  kitty::create_from_hex_string( function, "6996" );
  const auto cubes = mockturtle::exorcism( function );
  auto reconstructed = function.construct();
  kitty::create_from_cubes( reconstructed, cubes, true );
  return function == reconstructed;
}
