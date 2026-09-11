#include <mockturtle/mockturtle.hpp>

#ifndef MOCKTURTLE_PARENT_SENTINEL
#error Parent compile definitions were lost
#endif

int main()
{
  if ( fmt::format( FMT_STRING( "{}" ), 7 ) != "7" )
    return 1;
  mockturtle::aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  aig.create_po( aig.create_and( a, b ) );
  return aig.num_gates() == 1u && aig.num_pis() == 2u && aig.num_pos() == 1u ? 0 : 1;
}
