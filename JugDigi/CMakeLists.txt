################################################################################
# Package: JugDigi
################################################################################
gaudi_subdir(JugDigi v1r0)

find_package(EICD REQUIRED)
#set(PODIO $ENV{PODIO})
#set(CMAKE_MODULE_PATH CMAKE_MODULE_PATH PODIO)
#find_package(podio 0.11.01 REQUIRED)
#include_directories(${podio_INCLUDE_DIR})

find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)

# this declaration will not be needed in the future
gaudi_depends_on_subdirs(GaudiAlg GaudiKernel)
#gaudi_install_scripts()
#gaudi_install_python_modules()

#gaudi_add_library(JugBase
#  src/*.cpp 
#                  INCLUDE_DIRS EICD PODIO ROOT $ENV{HOME}/stow/podio/include
#                  LINK_LIBRARIES GaudiAlgLib GaudiKernel ROOT DD4hep::DDG4IO
#                  PUBLIC_HEADERS JugBase)
#target_link_libraries(JugBase 
#  podio::podioRootIO
#  )

gaudi_add_module(JugDigiPlugins
  src/components/CalorimeterHitDigi.cpp 
  src/components/CrystalEndcapsDigi.cpp 
  src/components/EMCalorimeterDigi.cpp 
  src/components/EcalTungstenSamplingDigi.cpp 
  src/components/ExampleCaloDigi.cpp 
  src/components/HadronicCalDigi.cpp 
  src/components/PhotoMultiplierDigi.cpp 
  src/components/SOIPIXTrackerDigi.cpp 
  src/components/UFSDTrackerDigi.cpp 
  src/components/ScFiBarrelCalDigi.cpp 
  LINK_LIBRARIES GaudiAlgLib GaudiKernel JugBase ROOT NPDet::DD4podIO EICD::eicd)
target_compile_options(JugDigiPlugins PRIVATE -Wno-suggest-override)

#gaudi_add_test(ProduceForReadTest
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/simple_producer.py)
gaudi_add_test(ExCaloDigi
               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
               FRAMEWORK tests/options/example_calodigi.py
               )
#gaudi_add_test(CheckReadCollectionSize
#               ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               COMMAND python FWCore/tests/scripts/check_coll_after_read.py
#               DEPENDS ReadTest)
