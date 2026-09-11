set(NAUTY_BUILD_SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/src")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/" DESTINATION "${NAUTY_BUILD_SOURCE_DIR}")
execute_process(COMMAND ${CMAKE_COMMAND} -E env "CC=${CMAKE_C_COMPILER}" "./configure" "--disable-popcnt"
  WORKING_DIRECTORY "${NAUTY_BUILD_SOURCE_DIR}"
  RESULT_VARIABLE nauty_configure_result OUTPUT_VARIABLE nauty_configure_output
  ERROR_VARIABLE nauty_configure_error)
if(NOT nauty_configure_result EQUAL 0)
  message(FATAL_ERROR "Failed to configure Nauty: ${nauty_configure_output}\n${nauty_configure_error}")
endif()
# Match nauty.a's object list; the other C files are standalone command-line tools.
add_library(mockturtle_nauty STATIC EXCLUDE_FROM_ALL
  ${NAUTY_BUILD_SOURCE_DIR}/nauty.c ${NAUTY_BUILD_SOURCE_DIR}/nautil.c
  ${NAUTY_BUILD_SOURCE_DIR}/nausparse.c ${NAUTY_BUILD_SOURCE_DIR}/naugraph.c
  ${NAUTY_BUILD_SOURCE_DIR}/schreier.c ${NAUTY_BUILD_SOURCE_DIR}/naurng.c)
set_target_properties(mockturtle_nauty PROPERTIES EXPORT_NAME _nauty POSITION_INDEPENDENT_CODE ON
  C_VISIBILITY_PRESET hidden EXCLUDE_FROM_ALL $<NOT:$<BOOL:${MOCKTURTLE_INSTALL}>>)
target_compile_options(mockturtle_nauty PRIVATE
  "$<$<AND:$<CONFIG:Release,RelWithDebInfo,MinSizeRel>,$<C_COMPILER_ID:GNU,Clang,AppleClang>>:-ffunction-sections;-fdata-sections>")
# Configure may select ABI macros needed by its generated headers.
file(READ "${NAUTY_BUILD_SOURCE_DIR}/makefile" nauty_makefile)
if(nauty_makefile MATCHES "(^|\n)CFLAGS=([^\n]*)")
  separate_arguments(nauty_flags UNIX_COMMAND "${CMAKE_MATCH_2}")
  foreach(flag IN LISTS nauty_flags)
    if(flag MATCHES "^(-m|-D|-U)")
      target_compile_options(mockturtle_nauty PRIVATE "${flag}")
    endif()
  endforeach()
endif()
if(MOCKTURTLE_ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT nauty_ipo_supported OUTPUT nauty_ipo_error LANGUAGES C)
  if(NOT nauty_ipo_supported)
    message(FATAL_ERROR "MOCKTURTLE_ENABLE_IPO is unsupported for Nauty: ${nauty_ipo_error}")
  endif()
  set_property(TARGET mockturtle_nauty PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
endif()
mockturtle_instrument(mockturtle_nauty)
