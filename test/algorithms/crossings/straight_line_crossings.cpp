#include <catch.hpp>

#include <mockturtle/algorithms/crossings/straight_line_crossings.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/cover.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rank_view.hpp>

#include <type_traits>

using namespace mockturtle;

TEMPLATE_TEST_CASE(
    "Empty network straight-line crossing number", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  rank_view<depth_view<TestType>> const rank_ntk{};

  straight_line_crossings_stats st{};

  straight_line_crossings( rank_ntk, {}, &st );

  CHECK( st.num_crossings == 0 );
}

TEMPLATE_TEST_CASE(
    "Simple network straight-line crossing number", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  TestType ntk{};

  auto const x1 = ntk.create_pi();
  auto const x2 = ntk.create_pi();
  auto const a1 = ntk.create_and( x1, x2 );
  ntk.create_po( a1 );

  rank_view<depth_view<TestType>> rank_ntk{ depth_view<TestType>{ ntk } };

  straight_line_crossings_stats st{};

  straight_line_crossings( rank_ntk, {}, &st );

  CHECK( st.num_crossings == 0 );

  // there should be no crossings in the network no matter the order of the inputs
  rank_ntk.swap( rank_ntk.get_node( x1 ), rank_ntk.get_node( x2 ) );

  CHECK( st.num_crossings == 0 );
}

TEMPLATE_TEST_CASE(
    "Three layer lattice straight-line crossing number", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  rank_view<depth_view<TestType>> rank_ntk{};

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

  straight_line_crossings_stats st{};

  straight_line_crossings( rank_ntk, {}, &st );

  CHECK( st.num_crossings == 4 );
}

TEMPLATE_TEST_CASE(
    "Straight-line crossing number with inputs in different levels", "[crossings]",
    aig_network )
{
  rank_view<depth_view<TestType>> rank_ntk{};

  SECTION( "Asymmetrically stacked nodes" )
  {
    // rank 1
    auto const x1 = rank_ntk.create_pi();
    auto const x2 = rank_ntk.create_pi();
    auto const x3 = rank_ntk.create_pi();

    // rank 2
    auto const x4 = rank_ntk.create_and( x2, x3 );

    // rank 3
    auto const x5 = rank_ntk.create_and( x1, x4 );

    // rank 4
    auto const x6 = rank_ntk.create_and( x2, x5 );

    // output
    rank_ntk.create_po( x6 );

    // validate rank positions
    REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x1 ) ) == 0 );
    REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x2 ) ) == 1 );
    REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x3 ) ) == 2 );
    REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x4 ) ) == 0 );
    REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x5 ) ) == 0 );
    REQUIRE( rank_ntk.rank_position( rank_ntk.get_node( x6 ) ) == 0 );

    straight_line_crossings_stats st{};

    straight_line_crossings( rank_ntk, {}, &st );

    CHECK( st.num_crossings == 1 );
  }
  SECTION( "Symmetrically stacked nodes" )
  {
    // rank 1
    auto const x1 = rank_ntk.create_pi();
    auto const x2 = rank_ntk.create_pi();

    // rank 2
    auto const x3 = rank_ntk.create_and( x1, x2 );
    auto const x4 = rank_ntk.create_and( !x1, !x2 );

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

    straight_line_crossings_stats st{};

    straight_line_crossings( rank_ntk, {}, &st );

    CHECK( st.num_crossings == 4 );
  }
}
