cmake_minimum_required(VERSION 3.2.0)

gaudi_subdir(Examples)

gaudi_add_module(Examples src/*.cpp
                 LINK_LIBRARIES GaudiAlg)
