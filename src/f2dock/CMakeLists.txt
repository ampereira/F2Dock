CMAKE_MINIMUM_REQUIRED (VERSION 2.6)


FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/f2dock/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
  ../../inc/f2dock
)


LINK_DIRECTORIES(
  $ENV{TACC_FFTW3_LIB}
)


SET(F2DOCK_LINK_LIBS 
  fast-clash
  fast-GB 
  fast-hydro
  fast-LJ 
  fast-PQ
  fast-resCont
  PG-range
  math
  utils
  vol
  fft-utils
  misc-ident
#  ${FFTW_LIB}
  fftw3
  c
  m
  pthread
)

if(HAVE_LIBMOL)
SET(F2DOCK_LINK_LIBS  ${F2DOCK_LINK_LIBS} 
  hbondFilter
  libmol
)
endif(HAVE_LIBMOL)


ADD_EXECUTABLE (F2Dock ${SOURCE_FILES} ${INCLUDE_FILES})


TARGET_LINK_LIBRARIES(F2Dock ${F2DOCK_LINK_LIBS})
