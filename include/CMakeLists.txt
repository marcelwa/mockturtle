add_library(mockturtle INTERFACE)
add_library(mockturtle::mockturtle ALIAS mockturtle)
target_compile_features(mockturtle INTERFACE cxx_std_17)
target_include_directories(mockturtle SYSTEM INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
# The headers instantiate large templates and fmt requires the compiler to read
# and emit UTF-8, so consumers need these too.
target_compile_options(mockturtle INTERFACE "$<$<CXX_COMPILER_ID:MSVC>:/EHsc;/bigobj;/utf-8>")
target_link_libraries(mockturtle INTERFACE mockturtle_dependencies mockturtle_abcsat mockturtle_abcesop)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
  target_link_libraries(mockturtle INTERFACE stdc++fs)
endif()

if(MOCKTURTLE_ENABLE_ABC)
  target_link_libraries(mockturtle INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/lib/abc_static/libabc.a> ${CMAKE_DL_LIBS})
  target_compile_definitions(mockturtle INTERFACE MOCKTURTLE_ENABLE_ABC)
endif()
