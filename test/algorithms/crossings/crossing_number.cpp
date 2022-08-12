#include <catch.hpp>

#include <mockturtle/algorithms/crossings/crossing_graph.hpp>
#include <mockturtle/algorithms/crossings/crossing_number.hpp>
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
    "Empty network crossing number", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  crossing_graph<TestType> const cg{ TestType{} };

  crossing_number_stats st{};

  crossing_number<TestType>( cg, {}, &st );

  CHECK( st.num_crossings == 0 );
}

TEMPLATE_TEST_CASE(
    "Simple network crossing number", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  TestType ntk{};

  auto const x1 = ntk.create_pi();
  auto const x2 = ntk.create_pi();
  auto const a1 = ntk.create_and( x1, x2 );
  ntk.create_po( a1 );

  crossing_graph<TestType> const cg{ ntk };

  crossing_number_stats st{};

  crossing_number<TestType>( cg, {}, &st );

  CHECK( st.num_crossings == 0 );
}

TEMPLATE_TEST_CASE(
    "Maj network crossing number", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  TestType ntk{};

  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const f1 = ntk.create_maj( a, b, c );
  auto const f2 = ntk.create_maj( d, e, f1 );
  auto const f3 = ntk.create_maj( a, d, f1 );
  auto const f4 = ntk.create_maj( f1, f2, f3 );
  ntk.create_po( f4 );

  crossing_graph<TestType> const cg{ ntk };

  crossing_number_stats st{};

  crossing_number<TestType>( cg, {}, &st );

  if constexpr ( std::is_same_v<TestType, aig_network> )
  {
    CHECK( st.num_crossings == 2 );
  }
  else if constexpr ( std::is_same_v<TestType, xag_network> )
  {
    CHECK( st.num_crossings == 1 );
  }
  else
  {
    CHECK( st.num_crossings == 0 );
  }
}

TEMPLATE_TEST_CASE(
    "Maj network crossing number with inverters and PO fan-outs", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  TestType ntk{};

  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();

  auto const f1 = ntk.create_maj( a, b, c );
  auto const f2 = ntk.create_maj( f1, c, d );
  ntk.create_po( f1 );
  ntk.create_po( !f1 );
  ntk.create_po( f2 );
  ntk.create_po( f2 );
  ntk.create_po( !f2 );

  crossing_graph<TestType> const cg{ ntk };

  crossing_number_stats st{};

  crossing_number<TestType>( cg, {}, &st );

  CHECK( st.num_crossings == 0 );
}
