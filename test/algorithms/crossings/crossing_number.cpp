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

TEMPLATE_TEST_CASE(
    "ZZZZ Crossing number with coordinate list embedding from rank_view", "[crossings]",
    aig_network /*, mig_network, xag_network, xmg_network, klut_network, cover_network */ )
{
  rank_view<depth_view<TestType>> ntk{};

  // rank 1
  auto const x1 = ntk.create_pi();
  auto const x2 = ntk.create_pi();
  auto const x3 = ntk.create_pi();

  // rank 2
  auto const x4 = ntk.create_and( x1, x2 );
  auto const x5 = ntk.create_and( x1, x3 );
  auto const x6 = ntk.create_and( x2, x3 );

  // rank 3
  auto const x7 = ntk.create_and( x4, x5 );
  auto const x8 = ntk.create_and( x4, x6 );
  auto const x9 = ntk.create_and( x5, x6 );

  // outputs
  ntk.create_po( x7 );
  ntk.create_po( x8 );
  ntk.create_po( x9 );

  // validate rank positions
  REQUIRE( ntk.rank_position( ntk.get_node( x1 ) ) == 0 );
  REQUIRE( ntk.rank_position( ntk.get_node( x2 ) ) == 1 );
  REQUIRE( ntk.rank_position( ntk.get_node( x3 ) ) == 2 );
  REQUIRE( ntk.rank_position( ntk.get_node( x4 ) ) == 0 );
  REQUIRE( ntk.rank_position( ntk.get_node( x5 ) ) == 1 );
  REQUIRE( ntk.rank_position( ntk.get_node( x6 ) ) == 2 );
  REQUIRE( ntk.rank_position( ntk.get_node( x7 ) ) == 0 );
  REQUIRE( ntk.rank_position( ntk.get_node( x8 ) ) == 1 );
  REQUIRE( ntk.rank_position( ntk.get_node( x9 ) ) == 2 );

  // construct crossing graph
  crossing_graph const cg{ ntk };

  crossing_number_params params{};
  crossing_number_stats st{};

  // specify coordinate list embedding
  params.initial_embedding = crossing_number_params::embedding_scheme::COORDINATE_LIST;

  crossing_number( cg, params, &st );

  CHECK( st.num_crossings == 4 );
}
