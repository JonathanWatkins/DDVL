PROJECT(test NONE)

# - Main CMake script for DDVL
#
# - Require lastest generatation of CMake
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# - Use the intel icl compiler
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	MESSAGE(STATUS ${CMAKE_SYSTEM_NAME}=LINUX)
	set(CMAKE_CXX_COMPILER "icc")
#	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2" )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	MESSAGE(STATUS ${CMAKE_SYSTEM_NAME}=NOTLINUX)
	set(CMAKE_CXX_COMPILER "icl")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi /O2 " )
  	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")	
	add_definitions(-D__WINDOWS__)
	add_definitions(-DTRILIBRARY)
	add_definitions(-DANSI_DECLARATORS)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")


PROJECT(DDVL CXX)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
# do nothing

ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
# set boost directories
	set(Boost_INCLUDE_DIR c:/boost_1_51_0)
	set(Boost_LIBRARY_DIR c:/boost_1_51_0/stage/lib)

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

# - Tell compiler where to find headers
include_directories(
	include
)

add_library(DDVLlib 
include/CSimulation.cpp
include/CParticle.cpp
include/CCell.cpp
include/CFileOutput.cpp 
include/CParallelEulerIntegrator.cpp
include/GeometryBase.cpp
#include/GeometryChannel.cpp
include/GeometryTube.cpp
include/rv_library.cpp
)


add_executable(DDVL DDVL.cpp)

# - Link DDVL to the DDVL library
target_link_libraries(DDVL DDVLlib)

# - Install DDVL to the bin subdirectory of CMAKE_INSTALL_PREFIX
#install(TARGETS DDVL DESTINATION bin)
