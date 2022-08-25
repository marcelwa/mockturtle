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