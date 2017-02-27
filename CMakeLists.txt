# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PRadSim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Tell CMake where to find the ROOT installation
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

#---------------------------------------------------------------------------
# Add thirdparty library and includes
if(NOT EXISTS $ENV{EVIODIR})
  include_directories(${PROJECT_SOURCE_DIR}/evio/include)
  link_directories(${PROJECT_SOURCE_DIR}/evio/lib64)
else()
  include_directories($ENV{EVIODIR}/include)
  link_directories($ENV{EVIODIR}/lib)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(PRadSim PRadSim.cc ${sources} ${headers})
target_link_libraries(PRadSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} evio evioxx)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B3. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(PRADSIM_SCRIPTS
  gui.mac run1.mac run2.mac vis.mac
  )

foreach(_script ${PRADSIM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS PRadSim DESTINATION bin)