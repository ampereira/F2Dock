#
# This macro is for setting up Boost for a target
#

macro(SetupBoost TargetName)
 #set(Boost_USE_STATIC_LIBS OFF)
 set(Boost_USE_MULTITHREADED ON)
 find_package(Boost 1.34.0 COMPONENTS thread date_time regex filesystem system)
 if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})
   set(LINK_LIBS 
     ${LINK_LIBS} ${Boost_LIBRARIES}
   ) 
   message("Boost includes: ${Boost_INCLUDE_DIRS}")
   message("Boost libraries: ${Boost_LIBRARIES}")
 else(Boost_FOUND)
   message(SEND_ERROR "If you're having trouble finding boost, set environment variables "
           "BOOST_INCLUDEDIR and BOOST_LIBRARYDIR to the appropriate paths")
 endif(Boost_FOUND)
 
 target_link_libraries(${TargetName}
   ${LINK_LIBS}
 )
endmacro(SetupBoost)
