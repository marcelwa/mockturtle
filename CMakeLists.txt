cmake_minimum_required(VERSION 3.21...4.0)
project(mockturtle LANGUAGES CXX)

# Options
option(MOCKTURTLE_BUILD_EXAMPLES "Build mockturtle examples" ${PROJECT_IS_TOP_LEVEL})
option(MOCKTURTLE_BUILD_TESTS "Build mockturtle tests" OFF)
option(MOCKTURTLE_BUILD_EXPERIMENTS "Build mockturtle experiments" OFF)
option(MOCKTURTLE_INSTALL "Install mockturtle's headers and CMake package" ${PROJECT_IS_TOP_LEVEL})
option(BILL_Z3 "Enable Z3 interface for bill library" OFF)
option(MOCKTURTLE_ENABLE_COVERAGE "Enable coverage reporting for gcc/clang" OFF)
option(MOCKTURTLE_ENABLE_MATPLOTLIB "Enable matplotlib library in experiments" OFF)
option(MOCKTURTLE_ENABLE_NAUTY "Enable the Nauty library for percy" OFF)
option(MOCKTURTLE_ENABLE_ABC "Enable linking ABC as a static library" OFF)
option(MOCKTURTLE_ENABLE_ASAN "Enable AddressSanitizer for mockturtle" OFF)

if((MOCKTURTLE_ENABLE_COVERAGE OR MOCKTURTLE_ENABLE_ASAN) AND NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  message(FATAL_ERROR "MOCKTURTLE_ENABLE_COVERAGE and MOCKTURTLE_ENABLE_ASAN require GCC or Clang")
endif()
if(WIN32)
  set(MOCKTURTLE_ENABLE_NAUTY OFF)
endif()

if(MOCKTURTLE_ENABLE_ABC AND NOT EXISTS ${PROJECT_SOURCE_DIR}/lib/abc_static/libabc.a)
  message(FATAL_ERROR "Cannot find libabc.a in lib/abc_static/. For more info, see https://github.com/lsils/abc-staticlib")
endif()

include(GNUInstallDirs)
include(cmake/MockturtleTargets.cmake)

# Vendored dependency header directories, shared by the build and the install rules.
set(MOCKTURTLE_HEADER_DIRECTORIES parallel-hashmap fmt kitty rang lorina json percy bill abcsat abcesop)

add_subdirectory(lib)
add_subdirectory(include)

if(MOCKTURTLE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(MOCKTURTLE_BUILD_TESTS)
  add_subdirectory(test)
endif()

if(MOCKTURTLE_BUILD_EXPERIMENTS)
  add_subdirectory(experiments)
endif()

if(MOCKTURTLE_INSTALL)
  include(cmake/MockturtleInstall.cmake)
endif()
