FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/fast-clash/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${HEADER_FILES})

INCLUDE_DIRECTORIES(

  ../../inc
  ../../inc/fast-clash
)

SET(POST_FILTERING_LINK_LIBS 

  fast-clash
  math
  utils
  pthread
)

ADD_EXECUTABLE (FilterOutput FilterOutput.cpp ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(FilterOutput ${POST_FILTERING_LINK_LIBS})
