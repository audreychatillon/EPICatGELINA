cmake_minimum_required (VERSION 3.0) 
# Setting the policy to match Cmake version
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

project(nptool_project CXX)                                                              
set(CMAKE_CXX_STANDARD 17)

add_compile_options("-Wno-return-type-c-linkage") 

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/bin )   
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib )   
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include ) 


### nptool ###
# find package 
find_package(nptool QUIET)
if(${nptool_FOUND})
  message(STATUS "nptool found: ${nptool_VERSION}" )
  # use this definition to isolate nptool dependent code
  add_compile_definitions(nptool_FOUND)
  # add include directory
  include_directories(${nptool_INCLUDE_DIRS})
  include_directories($ENV{NPTOOL_HOME}/$ENV{NPTOOL_ENV}/include)  
  # add include directory
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
else()
  message(STATUS "nptool not found: compiling with no nptool support" )
endif()

### ROOT ###
# find package 
find_package(ROOT QUIET)
if(${ROOT_FOUND})
  message(STATUS "ROOT found: ${ROOT_VERSION}" )
  # use this definition to isolate ROOT dependent code
  add_compile_definitions(ROOT_FOUND)
  # add usefull root fonction for dictionnary definition
  include(${ROOT_USE_FILE})
 else()
  message(STATUS "ROOT not found: compiling with no ROOT support" )
endif()

### GEANT4 ###
find_package(Geant4 QUIET)
if(${Geant4_FOUND})
  message(STATUS "Geant4 found: ${Geant4_VERSION}" )
  # use this definition to isolate Geant4 dependent code
  add_compile_definitions(Geant4_FOUND)
  # add usefull root fonction for dictionnary definition
  include(${Geant4_USE_FILE})
  # add include directory
  include_directories(${Geant4_INCLUDE_DIRS})
else()
  message(STATUS "Geant4 not found: compiling with no Geant4 support" )
endif()

if(EXISTS "duplicate.cmake")     
  include("duplicate.cmake")
endif()

add_subdirectory("analysis")
if(${nptool_Geant4_FOUND})
  add_subdirectory("simulation")
endif()
