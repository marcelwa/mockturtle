#include <catch.hpp>

#include <mockturtle/algorithms/crossings/crossing_graph.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/cover.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>

using namespace mockturtle;

template<typename Ntk, typename Fun>
void foreach_edge( Ntk const& ntk, Fun&& fun )
{
  ntk.foreach_node( [&ntk, &fun]( auto const& gate )
                    { ntk.foreach_fanin( gate, [&ntk, &fun, &gate]( auto const& fanin )
                                         {
                                           if (!ntk.is_constant(ntk.get_node(fanin)))
                                           {
                                             fun( gate, fanin );
                                           } } ); } );
}

template<typename Ntk>
void check_crossing_graph_edge_list( Ntk const& ntk, crossing_graph<Ntk> const& cg ) noexcept
{
  foreach_edge( ntk, [&ntk, edge_list = cg.get_edge_list(), M = cg.get_num_edges(), i = 0u]( auto const& gate, auto const& fanin ) mutable noexcept
                {
                    CHECK( edge_list[i] == static_cast<int>( ntk.node_to_index( gate ) ) );
                    CHECK( edge_list[i + M] == static_cast<int>( ntk.node_to_index( ntk.get_node( fanin ) ) ) );
                    ++i; } );
}

TEMPLATE_TEST_CASE(
    "Empty crossing_graph construction", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  crossing_graph<TestType> const cg{ TestType{} };

  CHECK( cg.get_num_vertices() == 0 );
  CHECK( cg.get_num_edges() == 0 );
}

TEMPLATE_TEST_CASE(
    "Simple crossing_graph construction", "[crossings]",
    aig_network, mig_network, xag_network, xmg_network, klut_network, cover_network )
{
  TestType ntk{};

  const auto x1 = ntk.create_pi();
  const auto x2 = ntk.create_pi();
  const auto a1 = ntk.create_and( x1, x2 );
  ntk.create_po( a1 );

  crossing_graph<TestType> const cg{ ntk };

  CHECK( cg.get_num_vertices() == 3 );
  CHECK( cg.get_num_edges() == 2 );

  check_crossing_graph_edge_list( ntk, cg );
}
