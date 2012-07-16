# Install script for directory: /Users/andre/F2Dock-refactored/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/Users/andre/F2Dock-refactored/src/fft-utils/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/math/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/misc-ident/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/utils/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/vol/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/fast-clash/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/fast-GB/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/fast-hydro/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/fast-LJ/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/fast-PQ/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/fast-resCont/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/PG-range/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/f2dock/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/GB-rerank/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/XmlRPC/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/F2DockServer/cmake_install.cmake")
  INCLUDE("/Users/andre/F2Dock-refactored/src/PostFiltering/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

