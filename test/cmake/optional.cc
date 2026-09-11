#ifdef CHECK_PLOT
#include <matplot/matplotlibcpp.h>
#include <filesystem>
#endif
#ifdef CHECK_Z3
#include <bill/sat/solver.hpp>
#endif
#ifdef CHECK_NAUTY
#include <percy/partial_dag.hpp>
#endif

int main()
{
#ifdef CHECK_NAUTY
  percy::partial_dag first( 2, 2 );
  first.set_vertex( 0, 0, 0 );
  first.set_vertex( 1, 0, 1 );
  const auto second = first;
  if ( !first.is_isomorphic( second ) )
    return 1;
#endif
#ifdef CHECK_Z3
  bill::solver<bill::solvers::z3> solver;
  const bill::lit_type literal( solver.add_variable(), bill::lit_type::polarities::positive );
  solver.add_clause( literal );
  if ( solver.solve() != bill::result::states::satisfiable )
    return 2;
  solver.add_clause( ~literal );
  if ( solver.solve() != bill::result::states::unsatisfiable )
    return 3;
#endif
#ifdef CHECK_PLOT
  matplotlibcpp::backend( "Agg" );
  if ( !matplotlibcpp::plot( std::vector<double>{ 0.0, 1.0 }, std::vector<double>{ 0.0, 1.0 } ) )
    return 4;
  matplotlibcpp::save( "optional-smoke.png" );
  if ( std::filesystem::file_size( "optional-smoke.png" ) == 0 )
    return 5;
#endif
  return 0;
}
