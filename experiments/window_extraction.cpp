//
// Created by marcel on 09.04.25.
//

#include <filesystem>
#include <string>

#include <kitty/hash.hpp>
#include <lorina/aiger.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/mockturtle.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/network_utils.hpp>
#include <mockturtle/utils/window_utils.hpp>
#include <mockturtle/views/color_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/topo_view.hpp>

#include <parallel_hashmap/phmap.h>

#include <experiments.hpp>

namespace mockturtle::hashing
{

template<typename VecTT>
struct hash
{
  std::size_t operator()( const VecTT& tts ) const
  {
    std::size_t seed = 0;
    for ( auto const& tt : tts )
    {
      kitty::hash_combine( seed, kitty::hash<kitty::dynamic_truth_table>()( tt ) );
    }
    return seed;
  }
};

} // namespace mockturtle::hashing

using namespace mockturtle;
namespace fs = std::filesystem;

int main()
{
  constexpr uint16_t MIN_CUT_SIZE = 6;
  constexpr uint16_t MAX_CUT_SIZE = 12;
  constexpr uint16_t MIN_NUM_LEVELS = 6;
  constexpr uint16_t MAX_NUM_LEVELS = 18;

  // Create directory structure if it doesn't exist
  fs::path windows_dir = "benchmarks/windows";
  fs::path aiger_dir = windows_dir / "aiger";
  fs::path dot_dir = windows_dir / "dot";

  if ( !fs::exists( windows_dir ) )
  {
    fs::create_directories( windows_dir );
  }
  if ( !fs::exists( aiger_dir ) )
  {
    fs::create_directories( aiger_dir );
  }
  if ( !fs::exists( dot_dir ) )
  {
    fs::create_directories( dot_dir );
  }

  // for ( auto const& benchmark : experiments::all_benchmarks() )
  // for ( auto const& benchmark : experiments::all_benchmarks( experiments::dec ) )
  for ( auto const& benchmark : experiments::all_benchmarks( experiments::iscas ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig{};
    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const fanout_view fanout_aig{ aig };
    const depth_view depth_aig{ fanout_aig };
    const color_view color_aig{ depth_aig };
    const topo_view topo_aig{ color_aig };

    // Extract benchmark name without path and extension
    const std::string benchmark_name = fs::path( benchmark ).stem().string();

    // Iterate over different cut sizes and level depths
    for ( auto cut_size = MIN_CUT_SIZE; cut_size <= MAX_CUT_SIZE; ++cut_size )
    {
      // truth table cache
      phmap::parallel_flat_hash_set<std::vector<kitty::dynamic_truth_table>, hashing::hash<std::vector<kitty::dynamic_truth_table>>> truth_table_cache{};

      for ( auto num_levels = MIN_NUM_LEVELS; num_levels <= MAX_NUM_LEVELS; ++num_levels )
      {
        fmt::print( "[i] processing with cut_size={}, num_levels={}\n", cut_size, num_levels );

        create_window_impl windowing( topo_aig );
        uint32_t window_index = 0;

        topo_aig.foreach_gate( [&topo_aig, &windowing, &window_index, &benchmark_name,
                                &aiger_dir, &dot_dir, cut_size, num_levels, &truth_table_cache]( const auto n ) {
          // window computation
          if ( const auto w = windowing.run( n, cut_size, num_levels ); w.has_value() )
          {
            // extract window as subnetwork
            aig_network window{};
            clone_subnetwork( topo_aig, w->inputs, w->outputs, w->nodes, window );

            // compute the window's truth tables
            const default_simulator<kitty::dynamic_truth_table> sim( window.num_pis() );
            const auto tts = simulate<kitty::dynamic_truth_table, aig_network>( window, sim );

            // if the truth tables are already in the cache, skip this window
            if ( truth_table_cache.contains( tts ) )
            {
              return;
            }

            // insert the truth tables into the cache
            truth_table_cache.insert( tts );

            // write the window to AIGER and DOT files
            const std::string window_filename = fmt::format( "{}_c_{}_l_{}_w_{}.aig",
                                                             benchmark_name, cut_size, num_levels, window_index );

            const std::string dot_filename = fmt::format( "{}_c_{}_l_{}_w_{}.dot",
                                                          benchmark_name, cut_size, num_levels, window_index );

            window_index++;

            const fs::path window_filepath = aiger_dir / window_filename;
            const fs::path dot_filepath = dot_dir / dot_filename;

            write_aiger( window, window_filepath.string() );
            write_dot( window, dot_filepath.string() );
          }
        } );
      }
    }
  }

  return 0;
}
