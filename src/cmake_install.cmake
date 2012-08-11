# Install script for directory: /h1/mmatos/F2Dock-refactored_02/src

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

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/fft-utils/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/math/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/misc-ident/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/utils/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/vol/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/fast-clash/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/fast-GB/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/fast-hydro/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/fast-LJ/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/fast-PQ/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/fast-resCont/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/PG-range/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/f2dock/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/GB-rerank/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/XmlRPC/cmake_install.cmake")
  INCLUDE("/h1/mmatos/F2Dock-refactored_02/src/F2DockServer/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

