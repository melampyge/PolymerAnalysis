cmake_minimum_required(VERSION 3.5.1)  
project(PolymerAnalysis)

#if(APPLE)
#  set(CMAKE_MACOSX_RPATH 1)
#  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
#  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#endif()

########## Version information ########## 

set (PolymerAnalysis_version_major 0)
set (PolymerAnalysis_version_minor 1)
set (PolymerAnalysis_version_patch 0)
set (CF_PACKAGE_VERSION ${PolymerAnalysis_version_major}.${PolymerAnalysis_version_minor}.${PolymerAnalysis_version_patch})

### Echo an empty line and program version

message("")
message(STATUS "configuring PolymerAnalysis version " ${PolymerAnalysis_version_major} . ${PolymerAnalysis_version_minor} . ${PolymerAnalysis_version_patch})

### Set cmake module directory

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules/" ${CMAKE_MODULE_PATH})
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

### Find necessary packages

find_package(HDF5 REQUIRED)
if(${HDF5_FOUND})
	message(STATUS "HDF5 is found with the properties: " ${HDF5_INCLUDE_DIRS})
	message(STATUS ${HDF5_LIBRARY_DIRS})
	message(STATUS ${HDF5_LIBRARIES})
	message(STATUS ${HDF5_DEFINITIONS})
else()
	message(FATAL_ERROR "HDF5 is not found!")
endif()

### Set compiler flags

add_definitions(-Wall -fno-inline -Wno-unused-variable -g -std=c++14 -O3 ${HDF5_DEFINITIONS})

### Set include directories

include_directories(include) 
include_directories(${HDF5_INCLUDE_DIRS})

### Set link directories

link_directories(${HDF5_LIBRARY_DIRS})

### Executable 

add_executable(exec src/analyser.cpp src/AnalyseNeighbours.cpp src/data_structures.cpp src/read_write.cpp)
target_link_libraries(exec ${HDF5_LIBRARIES})
#target_link_libraries(exec nsolver ${extralibs} ${PYTHON_LIBRARY} boostp boostn)

install(TARGETS exec DESTINATION ${CMAKE_SOURCE_DIR}/bin)

