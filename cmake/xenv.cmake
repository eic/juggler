###############################################################################
## Create and install minimal project xenv file
###############################################################################

configure_file("cmake/${PROJECT_NAME}.xenv.in"
  "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}.xenv"
  @ONLY)

install(FILES
  "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}.xenv"
  DESTINATION ${CMAKE_INSTALL_PREFIX} COMPONENT dev)
