cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Debug)


project(nsgaii)

find_package(MPI REQUIRED)

FIND_PACKAGE(Boost 1.4 COMPONENTS mpi serialization program_options REQUIRED)

INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} "main.cpp")


TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${MPI_LIBRARIES} ${Boost_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
