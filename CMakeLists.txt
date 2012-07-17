cmake_minimum_required (VERSION 2.6)
project (F2Dock)

# Global Variables
if(CMAKE_HOST_UNIX)
  # standard setup
  set (CMAKE_CXX_FLAGS "$ENV{CPPFLAGS} -I='/User/andre/F2Dock/inc/' -O2")
  set (CMAKE_LD_FLAGS "$ENV{LDFLAGS}")
  
  # Debugging and profiling friendly flags:
  #set (CMAKE_CXX_FLAGS "$ENV{CPPFLAGS} -pg")
  #set (CMAKE_LD_FLAGS "$ENV{LDFLAGS} -pg")
endif(CMAKE_HOST_UNIX)


# CMake Modules
SET(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_MODULE_PATH})

# include macro for setting up FFTW
include( SetupFFTW )

SetupFFTW()

#output 
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

IF(APPLE)
ADD_DEFINITIONS(-DMACOS_X)
ENDIF(APPLE)

add_subdirectory (src)
