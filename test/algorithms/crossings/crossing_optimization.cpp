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

  CHECK( st.num_crossings == 0 );
}