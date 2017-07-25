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

add_definitions(-Wall -fno-inline -Wno-unused-variable
-Wunused-but-set-variable -std=c++14 -O3 ${HDF5_DEFINITIONS})
if(UNIX)
		add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
endif()

### Set include directories

include_directories(include) 
include_directories(${HDF5_INCLUDE_DIRS})
if(UNIX)
		include_directories(/usr/local/gcc6/include/c++/6.3.0)
endif()

### Set link directories

link_directories(${HDF5_LIBRARY_DIRS})
if(UNIX)
		#add_library(stdc++ SHARED IMPORTED) 
		#set_target_properties(stdc++ PROPERTIES
  	#	IMPORTED_LOCATION "/usr/local/gcc6/lib64"
  	#	INTERFACE_INCLUDE_DIRECTORIES "/usr/local/gcc6/lib64"
		#)			
		link_directories(/usr/local/gcc6/lib64)
		link_directories(/usr/local/gcc6/lib/gcc/x86_64-redhat-linux/6.3.0)
endif()

### Executable 

add_executable(exec src/analyser.cpp src/AnalyseNeighbours.cpp src/data_structures.cpp src/read_write.cpp)
target_link_libraries(exec ${HDF5_LIBRARIES})
if(UNIX)
		target_link_libraries(exec stdc++ supc++ gcc_s)
endif()

install(TARGETS exec DESTINATION ${CMAKE_SOURCE_DIR}/bin)

