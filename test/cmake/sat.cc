#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/networks/aig.hpp>

static_assert( sizeof( pabc::ABC_PTRINT_T ) == sizeof( void* ) );
static_assert( sizeof( pabc::ABC_PTRUINT_T ) == sizeof( void* ) );

bool sat_check()
{
  mockturtle::equivalence_checking_params params;
  params.functional_reduction = false;
  mockturtle::aig_network lhs, rhs;
  const auto a = lhs.create_pi();
  const auto b = lhs.create_pi();
  lhs.create_po( lhs.create_and( a, b ) );
  const auto c = rhs.create_pi();
  const auto d = rhs.create_pi();
  rhs.create_po( rhs.create_and( d, c ) );
  const auto equal = mockturtle::equivalence_checking( *mockturtle::miter<mockturtle::aig_network>( lhs, rhs ), params );
  mockturtle::aig_network different;
  const auto e = different.create_pi();
  const auto f = different.create_pi();
  different.create_po( different.create_or( e, f ) );
  const auto unequal = mockturtle::equivalence_checking( *mockturtle::miter<mockturtle::aig_network>( lhs, different ), params );
  return equal && *equal && unequal && !*unequal;
}
