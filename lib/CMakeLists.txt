# The vendored headers are one target: mockturtle's public templates include
# them, and separate targets under generic names (`kitty`, `fmt`, ...) collide
# with a parent project's own targets and cannot be exported.
add_library(mockturtle_dependencies INTERFACE)
set_property(TARGET mockturtle_dependencies PROPERTY EXPORT_NAME dependencies)
foreach(directory IN LISTS MOCKTURTLE_HEADER_DIRECTORIES)
  target_include_directories(mockturtle_dependencies SYSTEM INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${directory}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mockturtle/dependencies/${directory}>)
endforeach()
target_compile_definitions(mockturtle_dependencies INTERFACE FMT_HEADER_ONLY ABC_NAMESPACE=pabc)

# ABC picks its integer widths from these. The previous hard-coded `LIN64` was
# wrong on any 32-bit target: `ABC_PTRINT_T` became wider than a pointer.
add_library(mockturtle_abc_platform INTERFACE)
set_property(TARGET mockturtle_abc_platform PROPERTY EXPORT_NAME abc_platform)
if(WIN32)
  target_compile_definitions(mockturtle_abc_platform INTERFACE
    $<IF:$<EQUAL:${CMAKE_SIZEOF_VOID_P},8>,NT64,NT>
    NOMINMAX WIN32_NO_DLL _CRT_SECURE_NO_WARNINGS)
elseif(CMAKE_SIZEOF_VOID_P EQUAL 8)
  target_compile_definitions(mockturtle_abc_platform INTERFACE LIN64)
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
  target_compile_definitions(mockturtle_abc_platform INTERFACE LIN)
else()
  message(FATAL_ERROR "Unsupported pointer size for ABC: ${CMAKE_SIZEOF_VOID_P}")
endif()
target_compile_definitions(mockturtle_abc_platform INTERFACE ABC_NO_USE_READLINE)
target_link_libraries(mockturtle_dependencies INTERFACE mockturtle_abc_platform)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
add_subdirectory(abcsat)
add_subdirectory(abcesop)

if(MOCKTURTLE_ENABLE_NAUTY)
  add_subdirectory(nauty)
  target_link_libraries(mockturtle_dependencies INTERFACE nauty)
  target_include_directories(mockturtle_dependencies SYSTEM INTERFACE
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mockturtle/dependencies/nauty>)
  target_compile_definitions(mockturtle_dependencies INTERFACE MOCKTURTLE_ENABLE_NAUTY)
else()
  target_compile_definitions(mockturtle_dependencies INTERFACE DISABLE_NAUTY)
endif()

if(BILL_Z3)
  set(BILL_Z3_INCLUDE_PATH "" CACHE PATH "Path to Z3 includes, e.g., z3++.h")
  set(BILL_Z3_LIBRARY_PATH "" CACHE PATH "Path to Z3 library, e.g., libz3.a")
  target_compile_definitions(mockturtle_dependencies INTERFACE BILL_HAS_Z3)
  if(NOT "${BILL_Z3_INCLUDE_PATH}" STREQUAL "")
    target_include_directories(mockturtle_dependencies SYSTEM INTERFACE ${BILL_Z3_INCLUDE_PATH})
  endif()
  if(NOT "${BILL_Z3_LIBRARY_PATH}" STREQUAL "")
    target_link_directories(mockturtle_dependencies INTERFACE ${BILL_Z3_LIBRARY_PATH})
  endif()
  target_link_libraries(mockturtle_dependencies INTERFACE $<IF:$<BOOL:${WIN32}>,libz3,z3>)
endif()

if(MOCKTURTLE_ENABLE_MATPLOTLIB)
  find_package(Python3 REQUIRED COMPONENTS Development NumPy)
  add_library(matplot INTERFACE)
  target_include_directories(matplot SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/matplot)
  target_link_libraries(matplot INTERFACE Python3::Python Python3::NumPy)
endif()
