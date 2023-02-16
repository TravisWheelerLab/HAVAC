cmake_minimum_required(VERSION 3.10)

# Set the name of your project
project(havac)

# Set the path to the source directory
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Set the build directory
set(BUILD_DIR "${SOURCE_DIR}/build")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 14)

# Set the C standard
set(CMAKE_C_STANDARD 11)

# Add directories to search for header files
include_directories(${SOURCE_DIR} ${XILINX_XRT}/include /usr/local/include)

# Add directories to search for libraries
link_directories(${XILINX_XRT}/lib /usr/local/lib)

# Add the source files for the executable target
target_sources(havac PRIVATE 
  ${SOURCE_DIR}/device/*.cpp
  ${SOURCE_DIR}/host/host/*.cpp
  ${SOURCE_DIR}/host/host/phmm/*.cpp
  ${SOURCE_DIR}/host/host/sequence/*.cpp
  ${SOURCE_DIR}/PhmmReprojection/*.cpp
)

# Create the executable target
add_executable(havac)

# Link against the required libraries
#xilinxopencl won't probably re required, but I'm keeping this here for the time being.
target_link_libraries(havac xilinxopencl xrt_coreutil pthread)
target_link_libraries(havac PUBLIC FastaVector P7HmmReader)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

# Set the output directory for the executable
set_target_properties(havac PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR})
