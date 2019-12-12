cmake_minimum_required(VERSION 3.16.0)

set(CMAKE_CXX_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/bin/intel64/icpc)
set(CMAKE_C_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/bin/intel64/icc)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_VERBOSE_MAKEFILE on)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -mkl=sequential")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mkl=sequential -ipo -inline-forceinline")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -restrict -mkl -ipo -inline-forceinline -no-prec-div -fp-model fast=2 -qopenmp")

set(MPI_CXX_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpiicpc)
set(MPI_CXX_HEADER_DIR /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include)

set(MPI_C_COMPILER /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/bin/mpiicc)
set(MPI_C_HEADER_DIR /opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include)

project(phynet LANGUAGES CXX VERSION 1.0)
#project(phynet LANGUAGES CXX CUDA VERSION 1.0)

include(CheckCCompilerFlag)
include(CheckIncludeFileCXX)

find_package(Git)
find_package(MKL REQUIRED)
#find_package(MPI REQUIRED)
#find_package(CUDA REQUIRED)
#find_package(OpenMP REQUIRED)
find_package(Boost COMPONENTS random REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE PATHS /usr/local/include/Eigen-3.3.7/cmake)
#find_package(ArrayFire REQUIRED PATHS /usr/local/include/arrayfire/share/ArrayFire/cmake)
	
add_definitions(-DEIGEN_USE_MKL_ALL)
set(EIGEN3_INCLUDE_DIRS /usr/local/include/Eigen-3.3.7)
set(ARRAY_FIRE_INCLUDE_DIRS /usr/local/include/arrayfire/include)

include_directories(../include ${EIGEN3_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}) 
file(GLOB SOURCES ../src/*.cxx)

add_executable(phynet.x ${SOURCES})
target_link_libraries(phynet.x Eigen3::Eigen ${Boost_RANDOM_LIBRARY} ${MKL_LIBRARIES})

#include_directories(../include ${EIGEN3_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS} ${MPI_CXX_HEADER_DIR} ${ARRAY_FIRE_INCLUDE_DIRS}) 
#file(GLOB SOURCES ../src/*.cxx)
#add_executable(phynetx ${SOURCES})
#target_link_libraries(phynet.x ArrayFire::af Eigen3::Eigen MPI::MPI_CXX ${Boost_RANDOM_LIBRARY} ${MKL_LIBRARIES} ${ArrayFire_LIBRARIES} ${CUDA_LIBRARIES})
