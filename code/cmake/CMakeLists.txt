cmake_minimum_required(VERSION 3.16.0)

#set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc)
#set(CMAKE_C_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/bin/intel64/icc)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_VERBOSE_MAKEFILE on)
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -restrict -mkl -ipo -inline-forceinline -no-prec-div -fp-model fast=2 -qopenmp")

#set(MPI_CXX_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpiicpc)
#set(MPI_CXX_HEADER_DIR /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include)

#set(MPI_C_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpiicc)
#set(MPI_C_HEADER_DIR /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include)

project(phynet LANGUAGES CXX VERSION 1.0)
#project(phynet LANGUAGES CXX CUDA VERSION 1.0)

include(CheckCCompilerFlag)
include(CheckIncludeFileCXX)

find_package(Git)
find_package(Boost REQUIRED COMPONENTS random)
find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})
add_definitions(-DHFD5PP_EIGEN)

#find_package(MPI REQUIRED)
#find_package(CUDA REQUIRED)
#find_package(OpenMP REQUIRED)
#find_package(ArrayFire REQUIRED PATHS /usr/local/include/arrayfire/share/ArrayFire/cmake)

find_package(Eigen3 3.3.7 REQUIRED NO_MODULE) 
include_directories("$ENV{EIGEN3_ROOT}")
include_directories("$ENV{EIGEN3_ROOT}/unsupported")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	find_package(MKL REQUIRED)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ipo -inline-forceinline")
	add_definitions(-DEIGEN_USE_MKL_ALL)
	include_directories(${MKL_INCLUDE_DIRS})
endif()
	
#set(EIGEN3_INCLUDE_DIR /usr/local/include/Eigen-3.3.7 /usr/local/include/Eigen-3.3.7/unsupported)
#set(ARRAY_FIRE_INCLUDE_DIRS /usr/local/include/arrayfire/include)

include_directories(../inc/phynet/ ../inc/datgen)  

add_executable(phynet.x 
	../../src/phynet/Activation.cxx  
	../../src/phynet/Dataset.cxx  
	../../src/phynet/Gradient.cxx  
	../../src/phynet/Layer.cxx  
	../../src/phynet/Loss.cxx  
	../../src/phynet/main.cxx  
	../../src/phynet/Model.cxx  
	../../src/phynet/Network.cxx  
	../../src/phynet/nixio.cxx  
	../../src/phynet/Optimizer.cxx  
	../../src/phynet/Parser.cxx  
	../../src/phynet/Prior.cxx  
	../../src/phynet/Topology.cxx)

add_executable(datgen.x 
	../../src/datgen/main.cxx
	../../src/datgen/Decomp.cxx
	../../src/datgen/operators.cxx
	../../src/datgen/pauli.cxx)

target_link_libraries(phynet.x ${Boost_LIBRARIES} ${MKL_LIBRARIES})
target_link_libraries(datgen.x ${Boost_LIBRARIES} ${MKL_LIBRARIES})

#include_directories(../include ${EIGEN3_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS} ${MPI_CXX_HEADER_DIR} ${ARRAY_FIRE_INCLUDE_DIRS}) 
#file(GLOB SOURCES ../src/*.cxx)
#add_executable(phynetx ${SOURCES})
#target_link_libraries(phynet.x ArrayFire::af Eigen3::Eigen MPI::MPI_CXX ${Boost_RANDOM_LIBRARY} ${MKL_LIBRARIES} ${ArrayFire_LIBRARIES} ${CUDA_LIBRARIES})
