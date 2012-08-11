#
# This macro is for setting up a sub-project to use CGAL
#

macro(SetupCGAL TargetName)
  if(NOT DISABLE_CGAL)
    find_package(CGAL)
    if(CGAL_FOUND)
      include(${CGAL_USE_FILE})
      # need the following flags in case CGAL has some special compiler needs for this compiler
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS_INIT}")
      set(LIBS ${LIBS} ${CGAL_LIBRARIES})
      add_definitions(-DUSING_CGAL)
      target_link_libraries(${TargetName} ${LIBS})
    else(CGAL_FOUND)
      message("${TargetName} is requesting CGAL but it isnt found on the system!")
    endif(CGAL_FOUND)
  endif(NOT DISABLE_CGAL)
endmacro(SetupCGAL)
