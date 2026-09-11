# The GIA wrapper declares the ABC API itself; only the namespace-enabled archive is needed.
# Build ABC as described by https://github.com/lsils/abc-staticlib.
find_library(MockturtleABC_LIBRARY NAMES abc libabc
  HINTS ${MockturtleABC_ROOT} PATH_SUFFIXES lib)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MockturtleABC REQUIRED_VARS MockturtleABC_LIBRARY)
if(MockturtleABC_FOUND AND NOT TARGET MockturtleABC::ABC)
  add_library(MockturtleABC::ABC UNKNOWN IMPORTED GLOBAL)
  set_target_properties(MockturtleABC::ABC PROPERTIES IMPORTED_LOCATION "${MockturtleABC_LIBRARY}")
endif()
mark_as_advanced(MockturtleABC_LIBRARY)
