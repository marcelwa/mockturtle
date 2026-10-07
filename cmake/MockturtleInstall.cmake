include(CMakePackageConfigHelpers)
set(MOCKTURTLE_PACKAGE_DIRECTORY ${CMAKE_INSTALL_LIBDIR}/cmake/mockturtle)

set(MOCKTURTLE_EXPORTED_TARGETS mockturtle mockturtle_dependencies mockturtle_abc_platform
                                mockturtle_abcsat mockturtle_abcesop)
if(MOCKTURTLE_ENABLE_NAUTY)
  list(APPEND MOCKTURTLE_EXPORTED_TARGETS nauty)
  install(DIRECTORY "${PROJECT_BINARY_DIR}/lib/nauty/src/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mockturtle/dependencies/nauty
    FILES_MATCHING PATTERN "*.h" PATTERN "COPYRIGHT")
endif()
install(TARGETS ${MOCKTURTLE_EXPORTED_TARGETS} EXPORT mockturtleTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY include/mockturtle DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp")
# Each vendored dependency keeps its own layout, so that its relative includes
# and its license file resolve exactly as they do in the source tree.
foreach(directory IN LISTS MOCKTURTLE_HEADER_DIRECTORIES)
  install(DIRECTORY "lib/${directory}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mockturtle/dependencies/${directory}"
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.inc"
    PATTERN "LICENSE*" PATTERN "COPYING*" PATTERN "COPYRIGHT*")
endforeach()
install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/mockturtle)

configure_package_config_file(cmake/mockturtleConfig.cmake.in
  ${PROJECT_BINARY_DIR}/mockturtleConfig.cmake
  INSTALL_DESTINATION ${MOCKTURTLE_PACKAGE_DIRECTORY})
install(FILES ${PROJECT_BINARY_DIR}/mockturtleConfig.cmake DESTINATION ${MOCKTURTLE_PACKAGE_DIRECTORY})
install(EXPORT mockturtleTargets NAMESPACE mockturtle:: DESTINATION ${MOCKTURTLE_PACKAGE_DIRECTORY})
