include_directories(catch2) # v2.2.1

file(GLOB_RECURSE FILENAMES *.cpp)
# package/ is a standalone CMake consumer driven by its own run.py, not a Catch2 suite.
list(FILTER FILENAMES EXCLUDE REGEX "/package/")

add_executable(run_tests ${FILENAMES})
target_link_libraries(run_tests PRIVATE mockturtle::mockturtle)
mockturtle_configure_executable(run_tests)
target_compile_definitions(run_tests PRIVATE BENCHMARKS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../experiments/benchmarks")
target_compile_definitions(run_tests PRIVATE CATCH_CONFIG_CONSOLE_WIDTH=300)
