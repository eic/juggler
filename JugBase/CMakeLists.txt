################################################################################
# Package: JugBase
################################################################################
gaudi_subdir(JugBase v1r0)

find_package(EICD)
find_package(podio)
find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)

# this declaration will not be needed in the future
gaudi_depends_on_subdirs(GaudiAlg GaudiKernel)
gaudi_install_scripts()
gaudi_install_python_modules()

gaudi_add_library(JugBase
		  src/*.cpp
                  INCLUDE_DIRS EICD PODIO ROOT
                  LINK_LIBRARIES GaudiAlgLib GaudiKernel
                  PUBLIC_HEADERS JugBase)
target_link_libraries(JugBase 
  podio::podioRootIO
  )

#gaudi_add_module(FWCorePlugins
#                 src/components/*.cpp
#                 LINK_LIBRARIES GaudiAlgLib GaudiKernel FWCore ROOT)
#
#gaudi_add_test(ProduceForReadTest
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/simple_producer.py)
#gaudi_add_test(ReadTest
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/simple_reader.py
#               DEPENDS ProduceForReadTest)
#gaudi_add_test(CheckReadCollectionSize
#               ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               COMMAND python FWCore/tests/scripts/check_coll_after_read.py
#               DEPENDS ReadTest)
