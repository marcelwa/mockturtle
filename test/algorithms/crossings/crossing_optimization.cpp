#include <catch.hpp>

#include <mockturtle/algorithms/crossings/crossing_optimization.hpp>
#include <mockturtle/algorithms/crossings/straight_line_crossings.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/cover.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/rank_view.hpp>

#include <type_traits>

using namespace mockturtle;

TEMPLATE_TEST_CASE(
    "Empty network crossing optimization", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  rank_view<depth_view<fanout_view<TestType>>> rank_ntk{};

  crossing_optimization( rank_ntk );

  CHECK( straight_line_crossings( rank_ntk ) == 0 );
}

TEMPLATE_TEST_CASE(
    "Simple network crossing optimization", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  rank_view<depth_view<fanout_view<TestType>>> rank_ntk{};

  auto const x1 = rank_ntk.create_pi();
  auto const x2 = rank_ntk.create_pi();
  auto const a1 = rank_ntk.create_and( x1, x2 );
  rank_ntk.create_po( a1 );

  crossing_optimization( rank_ntk );

  CHECK( straight_line_crossings( rank_ntk ) == 0 );
}

TEMPLATE_TEST_CASE(
    "Three layer lattice crossing optimization", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  rank_view<depth_view<fanout_view<TestType>>> rank_ntk{};

  // rank 1
  auto const x1 = rank_ntk.create_pi();
  auto const x2 = rank_ntk.create_pi();
  auto const x3 = rank_ntk.create_pi();

  // rank 2
  auto const x4 = rank_ntk.create_and( x1, x2 );
  auto const x5 = rank_ntk.create_and( x1, x3 );
  auto const x6 = rank_ntk.create_and( x2, x3 );

  // rank 3
  auto const x7 = rank_ntk.create_and( x4, x5 );
  auto const x8 = rank_ntk.create_and( x4, x6 );
  auto const x9 = rank_ntk.create_and( x5, x6 );

  // outputs
  rank_ntk.create_po( x7 );
  rank_ntk.create_po( x8 );
  rank_ntk.create_po( x9 );

  // validate rank positions
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x1 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x2 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x3 ) ) == 2 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x4 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x5 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x6 ) ) == 2 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x7 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x8 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x9 ) ) == 2 );

  // validate number of crossings prior to optimization
  REQUIRE( straight_line_crossings( rank_ntk ) == 4 );

  crossing_optimization( rank_ntk );

  // optimization should not be able to reduce the number of crossings in this instance
  CHECK( straight_line_crossings( rank_ntk ) == 4 );
}

TEMPLATE_TEST_CASE(
    "Crossing optimization with fanins in different levels", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  rank_view<depth_view<fanout_view<TestType>>> rank_ntk{};

  // rank 1
  auto const x1 = rank_ntk.create_pi();
  auto const x2 = rank_ntk.create_pi();

  // rank 2
  auto const x3 = rank_ntk.create_and( x1, x2 );
  auto const x4 = rank_ntk.create_or( x1, x2 );

  // rank 3
  auto const x5 = rank_ntk.create_and( x3, x2 );
  auto const x6 = rank_ntk.create_and( x4, x1 );

  // outputs
  rank_ntk.create_po( x5 );
  rank_ntk.create_po( x6 );

  // validate rank positions
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x1 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x2 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x3 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x4 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x5 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x6 ) ) == 1 );

  // 4 crossings prior to optimization
  REQUIRE( straight_line_crossings( rank_ntk ) == 4 );

  crossing_optimization( rank_ntk );

  // 1 crossing after
  CHECK( straight_line_crossings( rank_ntk ) == 1 );
}

TEMPLATE_TEST_CASE(
    "Crossing optimization with blocks", "[crossings]",
    klut_network, cover_network )
{
  rank_view<depth_view<fanout_view<TestType>>> rank_ntk{};

  // rank 1
  auto const x1 = rank_ntk.create_pi(); // 2
  auto const x2 = rank_ntk.create_pi(); // 3
  auto const x3 = rank_ntk.create_pi(); // 4

  // rank 2
  const auto n1 = rank_ntk.create_not( x1 );     // 5
  const auto n2 = rank_ntk.create_not( x2 );     // 6
  auto const x4 = rank_ntk.create_and( x2, x3 ); // 7

  // rank 3
  const auto n22 = rank_ntk.create_not( n2 );    // 8
  auto const x5 = rank_ntk.create_and( n1, x4 ); // 9

  // rank 4
  auto const x6 = rank_ntk.create_and( n22, x5 ); // 10

  // output
  rank_ntk.create_po( x6 );

  // validate rank positions
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x1 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x2 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x3 ) ) == 2 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( n1 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( n2 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x4 ) ) == 2 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( n22 ) ) == 0 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x5 ) ) == 1 );
  REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x6 ) ) == 0 );

  // 1 crossing prior to optimization
  REQUIRE( straight_line_crossings( rank_ntk ) == 1 );

  crossing_optimization( rank_ntk );

  // 0 crossings after
  CHECK( straight_line_crossings( rank_ntk ) == 0 );
}
