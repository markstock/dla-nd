#
# DLA-ND
# (c)2000-8,14-15 Mark J. Stock <markjstock@gmail.com>
#
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( dla-nd )

SET (CMAKE_CONFIGURATION_TYPES "Debug" "Release")
SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Select which configuration to build" )
SET (DO_STATIC FALSE CACHE BOOL "Link in libraries statically")
SET (CMAKE_VERBOSE_MAKEFILE on)
SET (CMAKE_EXPORT_COMPILE_COMMANDS on)

IF (CMAKE_COMPILER_IS_GNUCXX)
  # support Core2 at a minimum, haswell ideally
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c99 -march=core2 -mtune=haswell")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c99")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
ELSEIF (MSVC)
  SET (CMAKE_CXX_FLAGS "/Ox")
ENDIF ()

IF( WIN32 )

ELSEIF( APPLE )
  INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
  #FIND_LIBRARY(IOKIT_LIBRARY IOKit)
  #MARK_AS_ADVANCED ( IOKIT_LIBRARY)

  #INCLUDE_DIRECTORIES ( "/usr/local/include" )

  #SET( PLATFORM_LIBS ${FRAMEWORK_LIBS} ${IOKIT_LIBRARY} )
  SET( PLATFORM_LIBS m png z )

ELSEIF( UNIX )

  SET( PLATFORM_LIBS m png z )

ENDIF()

IF( DO_STATIC )
  SET (CPREPROCDEFS ${CPREPROCDEFS} -static)
ENDIF()

ADD_DEFINITIONS (${CPREPROCDEFS})
INCLUDE_DIRECTORIES( "src" )
SET(SOURCES "src/setup.c"
            "src/density.c"
            "src/inout.c"
            "src/ndtree.c"
            "src/main.c")

# create a binary for the 2D version
ADD_EXECUTABLE( "dla-2d" ${SOURCES} )
TARGET_LINK_LIBRARIES( "dla-2d" ${PLATFORM_LIBS} )
TARGET_COMPILE_DEFINITIONS( "dla-2d" PRIVATE -DDIM=2 -DNCHILD=4 )

# same for 3d and 4d
ADD_EXECUTABLE( "dla-3d" ${SOURCES} )
TARGET_LINK_LIBRARIES( "dla-3d" ${PLATFORM_LIBS} )
TARGET_COMPILE_DEFINITIONS( "dla-3d" PRIVATE -DDIM=3 -DNCHILD=8 )
ADD_EXECUTABLE( dla-4d ${SOURCES} )
TARGET_LINK_LIBRARIES( dla-4d ${PLATFORM_LIBS} )
TARGET_COMPILE_DEFINITIONS( dla-4d PRIVATE -DDIM=4 -DNCHILD=16 )

