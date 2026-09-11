# Consumed twice by run.py: once with mockturtle embedded through
# add_subdirectory, once against the installed package after its build tree and
# source tree are gone.
cmake_minimum_required(VERSION 3.21)
project(mockturtle_package_contract LANGUAGES CXX)

enable_testing()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS OFF CACHE BOOL "" FORCE)
add_compile_definitions(MOCKTURTLE_PARENT_SENTINEL)
get_directory_property(BEFORE_OPTIONS COMPILE_OPTIONS)
get_directory_property(BEFORE_DEFINITIONS COMPILE_DEFINITIONS)

# A parent project owns these names; mockturtle must not need them.
foreach(name kitty lorina fmt rang json percy bill parallel_hashmap libabcsat libabcesop)
  add_library(${name} INTERFACE)
endforeach()

if(MOCKTURTLE_SOURCE)
  add_subdirectory("${MOCKTURTLE_SOURCE}" mockturtle)
else()
  find_package(mockturtle CONFIG REQUIRED)
endif()

get_directory_property(AFTER_OPTIONS COMPILE_OPTIONS)
get_directory_property(AFTER_DEFINITIONS COMPILE_DEFINITIONS)
if(NOT CMAKE_CXX_STANDARD EQUAL 20 OR CMAKE_EXPORT_COMPILE_COMMANDS
   OR NOT BEFORE_OPTIONS STREQUAL AFTER_OPTIONS
   OR NOT BEFORE_DEFINITIONS STREQUAL AFTER_DEFINITIONS)
  message(FATAL_ERROR "mockturtle changed its parent's configuration")
endif()

add_executable(consumer consumer.cpp)
target_link_libraries(consumer PRIVATE mockturtle::mockturtle)
add_test(NAME consumer COMMAND consumer)
