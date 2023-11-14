cmake_minimum_required(VERSION 3.16)

project(HavacBenchmarks)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a build directory if it doesn't exist
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/build)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build)
endif()

# Set the build directory as the current binary directory
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)

add_executable(HavacBenchmarks HavacBenchmark.cpp)

target_link_libraries(HavacBenchmarks PRIVATE ${CMAKE_SOURCE_DIR}/../host/build/libhavac.a)
