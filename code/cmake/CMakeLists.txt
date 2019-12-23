cmake_minimum_required(VERSION 3.16.0)

project(phynet LANGUAGES CXX VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

#set(MPI_CXX_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpiicpc)
#set(MPI_CXX_HEADER_DIR /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include)

#set(MPI_C_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpiicc)
#set(MPI_C_HEADER_DIR /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include)

include(CheckCCompilerFlag)
include(CheckIncludeFileCXX)

find_package(Git)
find_package(Boost REQUIRED COMPONENTS random filesystem)

#find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ArrayFire NO_MODULE)

find_package(Eigen3 3.3.7 REQUIRED NO_MODULE) 
include_directories(../lib/eigen ../lib/eigen/unsupported)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	find_package(MKL REQUIRED)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ipo -inline-forceinline")
	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ipo -inline-forceinline -qoverride-limits -qopenmp")
	add_definitions(-DEIGEN_USE_MKL_ALL)
	include_directories(${MKL_INCLUDE_DIRS})
endif()
	
include_directories(../include)  

add_library(nixio STATIC ../../src/nixio/nixio.cxx)

add_executable(dump 
	../../src/dump/main.cxx
	../../src/dump/Reader.cxx
	../../src/gendat/Fields.cxx)

add_executable(phynet
	../../src/phynet/Activation.cxx  
	../../src/phynet/Dataset.cxx  
	../../src/phynet/Layer.cxx  
	../../src/phynet/Loss.cxx  
	../../src/phynet/main.cxx  
	../../src/phynet/Model.cxx  
	../../src/phynet/Network.cxx  
	../../src/gendat/Operators.cxx
	../../src/phynet/Optimizer.cxx  
	../../src/phynet/Parser.cxx  
	../../src/phynet/Topology.cxx)

add_executable(gendat 
	../../src/gendat/Generator.cxx
	../../src/gendat/Fields.cxx
	../../src/gendat/Instance.cxx
	../../src/gendat/main.cxx
	../../src/gendat/Operators.cxx
	../../src/gendat/Solver.cxx)

target_link_libraries(dump ${Boost_LIBRARIES})
target_link_libraries(gendat ${Boost_LIBRARIES} ${MKL_LIBRARIES})
target_link_libraries(phynet nixio ${Boost_LIBRARIES} ${MKL_LIBRARIES})





