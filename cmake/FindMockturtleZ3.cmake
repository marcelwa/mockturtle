find_path(MockturtleZ3_INCLUDE_DIR NAMES z3++.h
  HINTS ${MockturtleZ3_ROOT} ${BILL_Z3_INCLUDE_PATH} PATH_SUFFIXES include)
find_library(MockturtleZ3_LIBRARY NAMES z3 libz3
  HINTS ${MockturtleZ3_ROOT} ${BILL_Z3_LIBRARY_PATH} PATH_SUFFIXES lib bin)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MockturtleZ3 REQUIRED_VARS MockturtleZ3_INCLUDE_DIR MockturtleZ3_LIBRARY)
if(MockturtleZ3_FOUND AND NOT TARGET MockturtleZ3::Z3)
  add_library(MockturtleZ3::Z3 UNKNOWN IMPORTED GLOBAL)
  set_target_properties(MockturtleZ3::Z3 PROPERTIES
    IMPORTED_LOCATION "${MockturtleZ3_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${MockturtleZ3_INCLUDE_DIR}")
endif()
mark_as_advanced(MockturtleZ3_INCLUDE_DIR MockturtleZ3_LIBRARY)
