include(CMakeFindDependencyMacro)
find_dependency(podio REQUIRED)
find_dependency(edm4hep REQUIRED)

# - Include the targets file to create the imported targets that a client can
# link to (libraries) or execute (programs)
include("${CMAKE_CURRENT_LIST_DIR}/AlgorithmsTargets.cmake")
