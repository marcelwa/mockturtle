include(CMakePackageConfigHelpers)
set(mockturtle_package_directory ${CMAKE_INSTALL_LIBDIR}/cmake/mockturtle)
install(TARGETS mockturtle mockturtle_sat mockturtle_esop mockturtle_all
  mockturtle_dependencies mockturtle_abc_platform mockturtle_abcsat mockturtle_abcesop
  EXPORT mockturtleTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
if(MOCKTURTLE_ENABLE_NAUTY)
  install(TARGETS mockturtle_nauty EXPORT mockturtleTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY "${PROJECT_BINARY_DIR}/lib/nauty/src/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mockturtle/dependencies/nauty
    FILES_MATCHING PATTERN "*.h" PATTERN "COPYRIGHT")
endif()
install(DIRECTORY include/mockturtle DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.hpp")
# Keep each dependency's original layout, including relative includes and licenses.
set(mockturtle_header_directories parallel-hashmap fmt kitty rang lorina json percy bill abcsat abcesop)
if(MOCKTURTLE_ENABLE_MATPLOTLIB)
  list(APPEND mockturtle_header_directories matplot)
endif()
foreach(directory IN LISTS mockturtle_header_directories)
  install(DIRECTORY "lib/${directory}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mockturtle/dependencies/${directory}"
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.tpp" PATTERN "*.inc"
    PATTERN "LICENSE*" PATTERN "COPYING*" PATTERN "COPYRIGHT*")
endforeach()
install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/mockturtle)
configure_package_config_file(cmake/mockturtleConfig.cmake.in
  ${PROJECT_BINARY_DIR}/mockturtleConfig.cmake
  INSTALL_DESTINATION ${mockturtle_package_directory})
install(FILES ${PROJECT_BINARY_DIR}/mockturtleConfig.cmake
  cmake/FindMockturtleABC.cmake cmake/FindMockturtleZ3.cmake
  DESTINATION ${mockturtle_package_directory})
install(EXPORT mockturtleTargets NAMESPACE mockturtle:: DESTINATION ${mockturtle_package_directory})
