#
# This macro is for setting up a sub-project to use either Qt3 or Qt4
#

macro(SetupQt)
  # Find and setup Qt for this project.
  set(QT_REQUIRED TRUE)
  find_package(Qt)

  if(QT4_FOUND)
    find_package(Qt4 COMPONENTS QtCore QtGui QtXml QtOpenGL Qt3Support REQUIRED)
    set  (QT_USE_QTXML      TRUE)
    set  (QT_USE_QT3SUPPORT TRUE) 
    set  (QT_USE_QTOPENGL   TRUE)
    set  (QT_USE_QTCORE     TRUE)
    set  (QT_USE_QTGUI      TRUE)
    include(${QT_USE_FILE})
    set(QT3_FOUND FALSE)
    add_definitions(${QT_DEFINITIONS})
    add_definitions(-DQT_CLEAN_NAMESPACE)
    include_directories(${QT_INCLUDE_DIR})
    include_directories(${QT_QT_INCLUDE_DIR})
    include_directories(${QT_QTCORE_INCLUDE_DIR})
  elseif(QT_FOUND)
    set(QT_MT_REQUIRED TRUE)
    find_package(Qt3 REQUIRED)
    add_definitions(${QT_DEFINITIONS})
    include_directories(${QT_INCLUDE_DIR})
    set(QT3_FOUND TRUE)
  endif(QT4_FOUND)

  # Force using static version of Qt3 for now on Win32.
  # Ideally, we should check to see which version of the library is found
  # and set this define accordingly.  But we don't yet support DLLs for the
  # majority of VolRover, so it's not a big deal at the moment.
  if(WIN32)
    add_definitions(-DQT_NODLL)
  endif(WIN32)
endmacro(SetupQt)

macro(SetupQt3)
  message("SetupQt3 -- unused now")
endmacro(SetupQt3)
