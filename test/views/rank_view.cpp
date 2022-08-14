#include <catch.hpp>

#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/cover.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rank_view.hpp>

#include <memory>

using namespace mockturtle;

TEMPLATE_TEST_CASE( "Traits", "[rank_view]", aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  CHECK( is_network_type_v<TestType> );
  CHECK( !has_rank_position_v<TestType> );
  CHECK( !has_at_rank_position_v<TestType> );
  CHECK( !has_width_v<TestType> );
  CHECK( !has_foreach_node_in_rank_v<TestType> );

  using rank_ntk = rank_view<depth_view<TestType>>;

  CHECK( is_network_type_v<rank_ntk> );
  CHECK( has_rank_position_v<rank_ntk> );
  CHECK( has_at_rank_position_v<rank_ntk> );
  CHECK( has_width_v<rank_ntk> );
  CHECK( has_foreach_node_in_rank_v<rank_ntk> );

  using rank_rank_ntk = rank_view<rank_ntk>;

  CHECK( is_network_type_v<rank_rank_ntk> );
  CHECK( has_rank_position_v<rank_rank_ntk> );
  CHECK( has_at_rank_position_v<rank_rank_ntk> );
  CHECK( has_width_v<rank_rank_ntk> );
  CHECK( has_foreach_node_in_rank_v<rank_rank_ntk> );
}

TEMPLATE_TEST_CASE( "Compute ranks for a simple network", "[rank_view]", aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  TestType ntk{};
  auto const x1 = ntk.create_pi();
  auto const x2 = ntk.create_pi();
  auto const a1 = ntk.create_and( x1, x2 );
  ntk.create_po( a1 );

  depth_view depth_ntk{ ntk };

  rank_view rank_ntk{ depth_ntk };

  REQUIRE( rank_ntk.size() == ntk.size() );

  CHECK( rank_ntk.width() == 2u );

  auto const x1_lvl = rank_ntk.level( rank_ntk.get_node( x1 ) );
  auto const x2_lvl = rank_ntk.level( rank_ntk.get_node( x2 ) );
  auto const a1_lvl = rank_ntk.level( rank_ntk.get_node( a1 ) );

  CHECK( x1_lvl == 0u );
  CHECK( x2_lvl == 0u );
  CHECK( a1_lvl == 1u );

  auto const x1_pos = rank_ntk.rank_position( rank_ntk.get_node( x1 ) );
  auto const x2_pos = rank_ntk.rank_position( rank_ntk.get_node( x2 ) );
  auto const a1_pos = rank_ntk.rank_position( rank_ntk.get_node( a1 ) );

  CHECK( ( x1_pos == 0u || x1_pos == 1u ) );
  CHECK( ( x2_pos == 0u || x2_pos == 1u ) );
  CHECK( a1_pos == 0u );

  CHECK( rank_ntk.at_rank_position( x1_lvl, x1_pos ) == rank_ntk.get_node( x1 ) );
  CHECK( rank_ntk.at_rank_position( x2_lvl, x2_pos ) == rank_ntk.get_node( x2 ) );
  CHECK( rank_ntk.at_rank_position( a1_lvl, a1_pos ) == rank_ntk.get_node( a1 ) );

  rank_ntk.foreach_node_in_rank( 0ul, [&]( auto const& n, auto i )
                                 {
                                  switch (i)
                                  {
                                    case( 0ul ):
                                    {
                                        CHECK( n == rank_ntk.get_node( x1 ) );
                                        break;
                                    }
                                    case( 1ul ):
                                    {
                                        CHECK( n == rank_ntk.get_node( x2 ) );
                                        break;
                                    }
                                    default:
                                    {
                                      CHECK( false );
                                    }
                                   } } );

  rank_ntk.foreach_node_in_rank( 1ul, [&]( auto const& n )
                                 { CHECK( n == rank_ntk.get_node( a1 ) ); } );

  rank_ntk.swap( rank_ntk.get_node( x1 ), rank_ntk.get_node( x2 ) );

  CHECK( rank_ntk.at_rank_position( x1_lvl, x1_pos ) == rank_ntk.get_node( x2 ) );
  CHECK( rank_ntk.at_rank_position( x2_lvl, x2_pos ) == rank_ntk.get_node( x1 ) );

  rank_ntk.foreach_node_in_rank( 0ul, [&]( auto const& n, auto i )
                                 {
                                   switch (i)
                                   {
                                   case( 0ul ):
                                   {
                                     CHECK( n == rank_ntk.get_node( x2 ) );
                                     break;
                                   }
                                   case( 1ul ):
                                   {
                                     CHECK( n == rank_ntk.get_node( x1 ) );
                                     break;
                                   }
                                   default:
                                   {
                                     CHECK( false );
                                   }
                                   } } );
}

TEMPLATE_TEST_CASE( "compute ranks during node construction", "[rank_view]", aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  depth_view depth_ntk{ TestType{} };
  rank_view rank_ntk{ depth_ntk };

  auto const a = rank_ntk.create_pi();
  auto const b = rank_ntk.create_pi();
  auto const c = rank_ntk.create_pi();

  auto const a1 = rank_ntk.create_and( a, b );
  auto const a2 = rank_ntk.create_and( a1, c );
  rank_ntk.create_po( a2 );

  CHECK( rank_ntk.width() == 3u );

  rank_ntk.foreach_node_in_rank( 0ul, [&]( auto const& n, auto i )
                                 {
                                   switch (i)
                                   {
                                   case( 0ul ):
                                   {
                                     CHECK( n == rank_ntk.get_node( a ) );
                                     break;
                                   }
                                   case( 1ul ):
                                   {
                                     CHECK( n == rank_ntk.get_node( b ) );
                                     break;
                                   }
                                   case ( 2ul ):
                                   {
                                     CHECK( n == rank_ntk.get_node( c ) );
                                     break;
                                   }
                                   default:
                                   {
                                     CHECK( false );
                                   }
                                   } } );

  rank_ntk.foreach_node_in_rank( 1ul, [&]( auto const& n, auto i )
                                 {
                                   switch (i)
                                   {
                                   case( 0ul ):
                                   {
                                     CHECK( n == rank_ntk.get_node( a1 ) );
                                     break;
                                   }
                                   case( 1ul ):
                                   {
                                     CHECK( n == rank_ntk.get_node( a2 ) );
                                     break;
                                   }
                                   default:
                                   {
                                     CHECK( false );
                                   }
                                   } } );

  auto const a3 = rank_ntk.create_and( b, c );
  auto const a4 = rank_ntk.create_and( a, c );
  auto const o1 = rank_ntk.create_or( a, b );
  auto const o2 = rank_ntk.create_or( a3, a4 );
  rank_ntk.create_po( o1 );
  rank_ntk.create_po( o2 );

  CHECK( rank_ntk.width() == 4u );
}

TEMPLATE_TEST_CASE( "compute ranks during node construction after copy ctor", "[rank_view]", aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  rank_view<depth_view<TestType>> rank_ntk{};
  {
    auto tmp = std::make_unique<rank_view<depth_view<TestType>>>( rank_ntk );
    CHECK( rank_ntk.events().on_add.size() == 4u );

    rank_view cpy_rank{ *tmp }; // copy ctor
    CHECK( rank_ntk.events().on_add.size() == 6u );

    tmp.reset(); // don't access tmp anymore after this line!
    CHECK( rank_ntk.events().on_add.size() == 4u );

    auto const a = cpy_rank.create_pi();
    auto const b = cpy_rank.create_pi();
    auto const c = cpy_rank.create_pi();
    auto const t0 = cpy_rank.create_or( a, b );
    auto const t1 = cpy_rank.create_or( b, c );
    auto const t2 = cpy_rank.create_and( t0, t1 );
    auto const t3 = cpy_rank.create_or( b, t2 );
    cpy_rank.create_po( t3 );
    CHECK( cpy_rank.width() == 3u );

    auto const t4 = cpy_rank.create_and( a, c );
    auto const t5 = cpy_rank.create_or( a, c );
    auto const t6 = cpy_rank.create_and( t4, t5 );
    cpy_rank.create_po( t6 );
    CHECK( cpy_rank.width() == 4u );

    CHECK( rank_ntk.events().on_add.size() == 4u );
  }

  CHECK( rank_ntk.events().on_add.size() == 2u );
}

TEMPLATE_TEST_CASE( "compute ranks during node construction after copy assignment", "[rank_view]", aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  depth_view<TestType> xag{};
  rank_view dxag{ xag };
  {
    auto tmp = std::make_unique<rank_view<depth_view<TestType>>>( xag );
    dxag = *tmp; // copy assignment
    tmp.reset();
  }

  auto const a = dxag.create_pi();
  auto const b = dxag.create_pi();
  auto const c = dxag.create_pi();
  dxag.create_po( dxag.create_or( b, dxag.create_and( dxag.create_or( a, b ), dxag.create_or( b, c ) ) ) );

  CHECK( dxag.width() == 3u );
}
