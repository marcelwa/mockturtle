# Instrumentation requested through the options. Applied per target rather than
# to the directory, so that a project embedding mockturtle keeps its own flags.
function(mockturtle_instrument target)
  get_target_property(target_type ${target} TYPE)
  # A static archive links nothing itself; its consumers have to carry the
  # runtime flags instead, or the instrumented objects fail to link.
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

# Settings private to the vendored archives. They are compiled as
# position-independent code because consumers link them into shared libraries,
# such as Python extension modules.
function(mockturtle_configure_archive target)
  set_property(TARGET ${target} PROPERTY POSITION_INDEPENDENT_CODE ON)
  target_compile_features(${target} PRIVATE cxx_std_17)
  target_compile_options(${target} PRIVATE "$<$<CXX_COMPILER_ID:MSVC>:/EHsc;/bigobj;/utf-8>")
  mockturtle_instrument(${target})
endfunction()

# Warnings belong to mockturtle's own executables, never to its consumers.
function(mockturtle_configure_executable target)
  target_compile_options(${target} PRIVATE
    "$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wno-unknown-pragmas>"
    "$<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-gnu-anonymous-struct;-Wno-nested-anon-types>")
  mockturtle_instrument(${target})
endfunction()
