cmake_minimum_required(VERSION 3.10)

# Set the name of your project
project(havac)

# Set the path to the source directory
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Set the build directory
set(BUILD_DIR "${SOURCE_DIR}/build")

# set the shared library install directory
set(SHARED_LIB_INSTALL_DIR /usr/local/lib)

# set the static library install directory
set(STATIC_LIB_INSTALL_DIR ${BUILD_DIR})

# Create the build directory if it doesn't exist
if (NOT EXISTS ${BUILD_DIR})
    file(MAKE_DIRECTORY ${BUILD_DIR})
endif()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)


# Add directories to search for header files
include_directories(${SOURCE_DIR} /opt/xilinx/xrt/include /usr/local/include)

# Add directories to search for libraries
link_directories(/opt/xilinx/xrt/lib /usr/local/lib)

# Create the library target
add_library(havac)

# set the global header that should be installed
set_target_properties(havac PROPERTIES PUBLIC_HEADER "Havac.hpp;types/VerifiedHit.hpp")

# set where the library installs to
install(TARGETS havac
  ARCHIVE DESTINATION ${STATIC_LIB_INSTALL_DIR}
  LIBRARY DESTINATION ${SHARED_LIB_INSTALL_DIR}
  PUBLIC_HEADER DESTINATION /usr/local/include
)

# Add the source files for the library target
target_sources(havac PRIVATE 
  ${SOURCE_DIR}/host/Havac.cpp
  ${SOURCE_DIR}/host/HavacHwClient.cpp
  ${SOURCE_DIR}/host/HavacSSV.cpp
  ${SOURCE_DIR}/host/HitVerifier.cpp
  ${SOURCE_DIR}/host/phmm/PhmmPreprocessor.cpp
  ${SOURCE_DIR}/host/sequence/SequencePreprocessor.cpp
  ${SOURCE_DIR}/host/types/HardwareHitReport.cpp
  ${SOURCE_DIR}/host/types/VerifiedHit.cpp
  ${SOURCE_DIR}/PhmmReprojection/PhmmReprojection.cpp
)

# Link against the required libraries
#xilinxopencl won't probably be required, but I'm keeping this here for the time being.
target_link_libraries(havac PUBLIC xrt_coreutil uuid)
target_link_libraries(havac PUBLIC xilinxopencl pthread)
target_link_libraries(havac PUBLIC fastavector P7HmmReader)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
