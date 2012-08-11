CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/hbondFilter/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
  ../../inc/hbondFilter
)

LINK_DIRECTORIES(
  $ENV{TACC_FFTW3_LIB}
)

ADD_LIBRARY (hbondFilter STATIC ${SOURCE_FILES} ${INCLUDE_FILES})

SET(HBOND_LINK_LIBS 
  hbondFilter
  libmol
  fft-utils
  fftw3
  pthread
)

#ADD_EXECUTABLE (testHbondFilter testHbondFilter.cpp)

#TARGET_LINK_LIBRARIES(testHbondFilter ${HBOND_LINK_LIBS})
