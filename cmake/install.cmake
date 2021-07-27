###############################################################################
## Create and install minimal project xenv file
###############################################################################

configure_file("cmake/${PROJECT_NAME}.xenv.in"
  "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}.xenv"
  @ONLY)

install(FILES
  "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}.xenv"
  DESTINATION ${CMAKE_INSTALL_PREFIX} COMPONENT dev)

###############################################################################
## Create and install juggler executable, both for build tree
## and for install prefix
###############################################################################
## build-tree
set (JUGGLER_LIB_PATH ${PROJECT_BINARY_DIR})
set (JUGGLER_PYTHONPATH ${PROJECT_BINARY_DIR}/python)
configure_file("cmake/juggler.in"
  "${PROJECT_BINARY_DIR}/juggler"
  @ONLY)

## install tree
set (JUGGLER_LIB_PATH ${CMAKE_INSTALL_PREFIX}/lib)
set (JUGGLER_PYTHONPATH ${CMAKE_INSTALL_PREFIX}/python)
configure_file("cmake/juggler.in"
  "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/juggler"
  @ONLY)
install(PROGRAMS
  "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/juggler"
  DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime)


