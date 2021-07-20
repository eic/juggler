################################################################################
# Package: JugBase
################################################################################
gaudi_subdir(JugBase v1r0)

find_package(EICD)

find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )
find_package(ROOT COMPONENTS RIO Tree Core REQUIRED)
find_package(DD4hep COMPONENTS DDG4 DDG4IO DDRec REQUIRED)
find_library(genfit2 genfit2 /usr/local/lib REQUIRED)


# this declaration will not be needed in the future
gaudi_depends_on_subdirs(GaudiAlg GaudiKernel)
#gaudi_install_scripts()
#gaudi_install_python_modules()

gaudi_add_library(JugBase
  src/*.cpp 
  src/Utilities/*.cpp
  src/Plugins/BFieldOptions.cpp  
  #src/Plugins/BFieldScalor.cpp  
  src/Plugins/BFieldUtils.cpp
                  INCLUDE_DIRS EICD PODIO ROOT $ENV{HOME}/stow/podio/include
                  LINK_LIBRARIES GaudiAlgLib GaudiKernel ROOT DD4hep::DDG4IO
                  PUBLIC_HEADERS JugBase )
target_link_libraries(JugBase 
  podio::podioRootIO
  ActsCore ActsPluginDD4hep
  )
target_compile_options(JugBase PRIVATE -Wno-suggest-override)

gaudi_add_module(JugBasePlugins
                 src/components/*.cpp 
                 src/components/ReadTestConsumer.cxx
                 LINK_LIBRARIES GaudiAlgLib GaudiKernel JugBase ROOT DD4hep::DDRec NPDet::DD4podIO  ActsCore ActsPluginDD4hep EICD::eicd genfit2)
target_compile_options(JugBasePlugins PRIVATE -Wno-suggest-override)

#gaudi_add_test(ProduceForReadTest
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               FRAMEWORK tests/options/simple_producer.py)
gaudi_add_test(ReadTest
               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
               FRAMEWORK tests/options/simple_reader.py)
gaudi_add_test(ReadGeoTest
               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
               FRAMEWORK tests/options/reader_with_geosvc.py)
#gaudi_add_test(CheckReadCollectionSize
#               ENVIRONMENT PYTHONPATH+=${PODIO_PYTHON_DIR}
#               WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#               COMMAND python JugBase/tests/scripts/check_coll_after_read.py
#               DEPENDS ReadTest)
