function(mockturtle_instrument target)
  get_target_property(target_type ${target} TYPE)
  if(target_type STREQUAL "STATIC_LIBRARY")
    set(link_scope INTERFACE)
  else()
    set(link_scope PRIVATE)
  endif()
  if(MOCKTURTLE_ENABLE_COVERAGE)
    target_compile_options(${target} PRIVATE -O0 -g --coverage)
    target_link_options(${target} ${link_scope} --coverage)
  endif()
  if(MOCKTURTLE_ENABLE_ASAN)
    target_compile_options(${target} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
    target_compile_definitions(${target} PRIVATE ADDRESS_SANITIZER)
    target_link_options(${target} ${link_scope} -fsanitize=address)
  endif()
endfunction()

# Settings private to the embedded compiled libraries.
function(mockturtle_configure_archive target)
  set_target_properties(${target} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    EXCLUDE_FROM_ALL $<NOT:$<BOOL:${MOCKTURTLE_INSTALL}>>)
  target_compile_features(${target} PRIVATE cxx_std_17)
  if(MOCKTURTLE_ENABLE_IPO)
    set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
  endif()
  target_compile_options(${target} PRIVATE
    "$<$<AND:$<CONFIG:Release,RelWithDebInfo,MinSizeRel>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-ffunction-sections;-fdata-sections>"
    "$<$<AND:$<CONFIG:Release,RelWithDebInfo,MinSizeRel>,$<CXX_COMPILER_ID:MSVC>>:/Gy;/Gw>"
    "$<$<CXX_COMPILER_ID:MSVC>:/EHsc;/bigobj;/utf-8>")
  mockturtle_instrument(${target})
endfunction()

# Warnings belong to mockturtle's own executables, never its consumers.
function(mockturtle_configure_executable target)
  target_compile_options(${target} PRIVATE
    "$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra>"
    "$<$<CXX_COMPILER_ID:MSVC>:/EHsc;/bigobj;/utf-8>")
  mockturtle_instrument(${target})
endfunction()

function(mockturtle_find_optional_dependencies)
  list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
  if(MOCKTURTLE_ENABLE_ABC)
    find_package(MockturtleABC REQUIRED)
  endif()
  if(BILL_Z3)
    find_package(MockturtleZ3 REQUIRED)
  endif()
  if(MOCKTURTLE_ENABLE_MATPLOTLIB)
    find_package(Python3 GLOBAL REQUIRED COMPONENTS Development NumPy)
  endif()
endfunction()
