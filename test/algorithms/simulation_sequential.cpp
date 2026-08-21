#include <catch.hpp>

#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/algorithms/simulation_sequential.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/sequential.hpp>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include <cstdint>
#include <string>
#include <vector>

using namespace mockturtle;

namespace
{

/*! \brief Builds a Fibonacci LFSR with taps on the top two bits.
 *
 * It has no primary inputs at all, so it runs off its reset state alone -- which
 * makes it a direct test of whether the reset values are honoured.
 */
sequential<aig_network> lfsr( uint32_t width, uint32_t seed )
{
  sequential<aig_network> aig;

  std::vector<aig_network::signal> state( width );
  for ( auto i = 0u; i < width; ++i )
  {
    state[i] = aig.create_ro();
  }

  auto const feedback = aig.create_xor( state[width - 1], state[width - 2] );

  /* primary outputs are created before register inputs: both are combinational
     outputs of the same network, sliced by position */
  aig.create_po( state[width - 1] );

  aig.create_ri( feedback );
  for ( auto i = 0u; i + 1 < width; ++i )
  {
    aig.create_ri( state[i] );
  }

  for ( auto i = 0u; i < width; ++i )
  {
    mockturtle::register_t reg;
    reg.init = ( ( seed >> i ) & 1 ) ? register_init::one : register_init::zero;
    aig.set_register( i, reg );
  }

  return aig;
}

/*! \brief Collects the single primary output of every cycle into a bit string. */
std::string trace_of( std::vector<std::vector<bool>> const& trace )
{
  std::string bits;
  for ( auto const& outputs : trace )
  {
    bits += outputs[0] ? '1' : '0';
  }
  return bits;
}

} /* namespace */

TEST_CASE( "simulate an LFSR from its reset state", "[simulation_sequential]" )
{
  auto const aig = lfsr( 4, 1 );

  auto const trace = simulate_sequential<bool>( aig, 15, default_simulator<bool>( std::vector<bool>{} ) );

  CHECK( trace.size() == 15 );

  /* a maximal-length sequence: 15 states before it comes back around */
  CHECK( trace_of( trace ) == "000100110101111" );

  /* and it does come back around -- cycle 15 repeats cycle 0 */
  auto const two_periods = simulate_sequential<bool>( aig, 30, default_simulator<bool>( std::vector<bool>{} ) );
  CHECK( trace_of( two_periods ).substr( 0, 15 ) == trace_of( two_periods ).substr( 15 ) );
}

TEST_CASE( "a different seed shifts the same sequence", "[simulation_sequential]" )
{
  /* seeding with the second state of the first LFSR must produce the same
     sequence one step ahead, which is only true if the reset values are used */
  auto const from_one = simulate_sequential<bool>( lfsr( 4, 1 ), 15, default_simulator<bool>( std::vector<bool>{} ) );
  auto const from_two = simulate_sequential<bool>( lfsr( 4, 2 ), 15, default_simulator<bool>( std::vector<bool>{} ) );

  CHECK( trace_of( from_one ).substr( 1 ) == trace_of( from_two ).substr( 0, 14 ) );
}

TEST_CASE( "a register with no reset value follows the parameter", "[simulation_sequential]" )
{
  /* a single register that simply holds whatever it was reset to */
  sequential<aig_network> aig;
  auto const state = aig.create_ro();
  aig.create_po( state );
  aig.create_ri( state );

  mockturtle::register_t reg;
  reg.init = register_init::unknown;
  aig.set_register( 0, reg );

  simulate_sequential_params ps;

  ps.undefined_reset_value = false;
  CHECK( trace_of( simulate_sequential<bool>( aig, 3, default_simulator<bool>( std::vector<bool>{} ), ps ) ) == "000" );

  ps.undefined_reset_value = true;
  CHECK( trace_of( simulate_sequential<bool>( aig, 3, default_simulator<bool>( std::vector<bool>{} ), ps ) ) == "111" );
}

TEST_CASE( "simulate a shift register with a per-cycle stimulus", "[simulation_sequential]" )
{
  /* three registers in a chain: whatever is put in appears at the output three
     cycles later */
  sequential<aig_network> aig;

  auto const in = aig.create_pi();
  auto const a = aig.create_ro();
  auto const b = aig.create_ro();
  auto const c = aig.create_ro();

  aig.create_po( c );

  aig.create_ri( in );
  aig.create_ri( a );
  aig.create_ri( b );

  for ( auto i = 0u; i < 3u; ++i )
  {
    mockturtle::register_t reg;
    reg.init = register_init::zero;
    aig.set_register( i, reg );
  }

  /* a single 1 on the input, then silence */
  stimulus_simulator sim( { { true }, { false } } );

  CHECK( trace_of( simulate_sequential<bool>( aig, 6, sim ) ) == "000100" );
}

TEST_CASE( "a stimulus shorter than the run holds its last assignment", "[simulation_sequential]" )
{
  sequential<aig_network> aig;

  auto const in = aig.create_pi();
  auto const state = aig.create_ro();

  aig.create_po( state );
  aig.create_ri( in );

  mockturtle::register_t reg;
  reg.init = register_init::zero;
  aig.set_register( 0, reg );

  /* One assignment for a four-cycle run: the input stays high after cycle 0.
     Spelled through a named vector rather than as `sim( { { true } } )`, which
     GCC 12 and older cannot tell apart from a copy construction -- the same
     reason `default_simulator<bool>` is spelled with an explicit `std::vector<bool>{}`
     throughout this file. */
  std::vector<std::vector<bool>> const stimulus{ { true } };
  stimulus_simulator sim( stimulus );

  CHECK( trace_of( simulate_sequential<bool>( aig, 4, sim ) ) == "0111" );
}

TEST_CASE( "simulate a sequential network with truth tables", "[simulation_sequential]" )
{
  /* a register holding the AND of the two primary inputs: the output is constant
     0 in the first cycle and the AND from the second one on */
  sequential<aig_network> aig;

  auto const x0 = aig.create_pi();
  auto const x1 = aig.create_pi();
  auto const state = aig.create_ro();

  aig.create_po( state );
  aig.create_ri( aig.create_and( x0, x1 ) );

  mockturtle::register_t reg;
  reg.init = register_init::zero;
  aig.set_register( 0, reg );

  auto const trace = simulate_sequential<kitty::dynamic_truth_table>(
      aig, 3, default_simulator<kitty::dynamic_truth_table>( 2 ) );

  kitty::dynamic_truth_table expected( 2 );
  kitty::create_from_hex_string( expected, "8" );

  CHECK( kitty::is_const0( trace[0][0] ) );
  CHECK( trace[1][0] == expected );
  CHECK( trace[2][0] == expected );
}

TEST_CASE( "simulating no cycles yields no values", "[simulation_sequential]" )
{
  CHECK( simulate_sequential<bool>( lfsr( 4, 1 ), 0, default_simulator<bool>( std::vector<bool>{} ) ).empty() );
}
