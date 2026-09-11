add_library(mockturtle_abcsat STATIC EXCLUDE_FROM_ALL
  AbcGlucose.cpp Glucose.cpp SimpSolver.cpp satSolver.cpp satStore.cpp)
set_property(TARGET mockturtle_abcsat PROPERTY EXPORT_NAME _abcsat)
target_include_directories(mockturtle_abcsat PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(mockturtle_abcsat PRIVATE ABC_NAMESPACE=pabc)
target_link_libraries(mockturtle_abcsat PUBLIC mockturtle_abc_platform Threads::Threads)
mockturtle_configure_archive(mockturtle_abcsat)
