#
# This macro is for setting up a sub-project to use the Gnu Scientific Library (GSL)
#

macro(SetupGSL TargetName)
  find_package(GSL)
  if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_GSL_CXX_FLAGS}")
    set(LIBS ${LIBS} ${GSL_LIBRARIES})
  else(GSL_FOUND)
    message(SEND_ERROR "${TargetName} requires the Gnu Scientific Library (GSL)!")
  endif(GSL_FOUND)
  target_link_libraries(${TargetName} ${LIBS})
endmacro(SetupGSL)
