#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/exorcism.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/networks/aig.hpp>

#include <kitty/constructors.hpp>
#include <parallel_hashmap/phmap.h>

#ifndef MOCKTURTLE_PARENT_SENTINEL
#error Parent compile definitions were dropped
#endif

// ABC derives its integer widths from the platform defines the package carries.
static_assert( sizeof( pabc::ABC_PTRINT_T ) == sizeof( void* ) );
static_assert( sizeof( abc::exorcism::ABC_PTRINT_T ) == sizeof( void* ) );

int main()
{
  // The ESOP backend.
  kitty::dynamic_truth_table function( 4u );
  kitty::create_from_hex_string( function, "6996" );
  auto reconstructed = function.construct();
  kitty::create_from_cubes( reconstructed, mockturtle::exorcism( function ), true );
  if ( function != reconstructed )
    return 1;

  // The SAT backend.
  mockturtle::aig_network lhs, rhs;
  const auto a = lhs.create_pi();
  const auto b = lhs.create_pi();
  lhs.create_po( lhs.create_and( a, b ) );
  const auto c = rhs.create_pi();
  const auto d = rhs.create_pi();
  rhs.create_po( rhs.create_and( d, c ) );

  mockturtle::equivalence_checking_params ps;
  ps.functional_reduction = false;
  const auto equivalent = mockturtle::equivalence_checking( *mockturtle::miter<mockturtle::aig_network>( lhs, rhs ), ps );
  if ( !equivalent || !*equivalent )
    return 2;

  phmap::flat_hash_map<int, int> map{ { 1, 2 } };
  return map.at( 1 ) == 2 ? 0 : 3;
}
