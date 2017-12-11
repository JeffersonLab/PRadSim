#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PRadSim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Tell CMake where to find the ROOT installation
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathMore Foam)

#---------------------------------------------------------------------------
# Locate pradanalyzer library and headers
#
if(DEFINED ENV{PRAD_PATH})
  set(PRADANA_INCLUDE_DIR $ENV{PRAD_PATH}/include CACHE PATH "Path to prad analyzer headers")
  set(PRADANA_LIBRARY_DIR $ENV{PRAD_PATH}/lib CACHE PATH "Path to prad analyzer libraries")
else()
  set(PRADANA_INCLUDE_DIR PRADANA_INCLUDE_DIR_NOT_FOUND CACHE PATH "Path to prad analyzer headers")
  set(PRADANA_LIBRARY_DIR PRADANA_INCLUDE_DIR_NOT_FOUND CACHE PATH "Path to prad analyzer libraries")
endif()
set(PRADANA_LIBRARIES prana prconf cana)

#----------------------------------------------------------------------------
# Setup include and lib directory for this project
#
if(NOT ${PRADANA_INCLUDE_DIR} STREQUAL "PRADANA_INCLUDE_DIR_NOT_FOUND" OR ${PRADANA_LIBRARY_DIR} STREQUAL "PRADANA_LIBRARY_DIR_NOT_FOUND")
  include_directories(${PRADANA_INCLUDE_DIR})
  link_directories(${PRADANA_LIBRARY_DIR})
else()
  message(FATAL_ERROR "PRad Analyzer headers and libraries not found")
endif()
include(${Geant4_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB PRADSIM_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB PRADSIM_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries
#
add_executable(PRadSim PRadSim.cc ${PRADSIM_SOURCES} ${PRADSIM_HEADERS})
target_link_libraries(PRadSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${PRADANA_LIBRARIES})

#----------------------------------------------------------------------------
# Add subdirectories to build
#
add_subdirectory(digitization)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory
# This is so that we can run the executable directly because it relies on these scripts being in the current working directory
#
set(PRADSIM_RESOURCE_FILES init_vis.mac vis.mac gui.mac run.mac)

add_custom_target(copy_files ALL COMMENT "Copy files")
foreach(_file ${PRADSIM_RESOURCE_FILES})
  add_custom_command(TARGET copy_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${_file} ${PROJECT_BINARY_DIR}/${_file}
  )
endforeach()

add_custom_target(link_prana ALL COMMENT "Link PRadAnalyzer files" DEPENDS PRadSim)
foreach(_lib ${PRADANA_LIBRARIES})
  add_custom_command(TARGET link_prana PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${PRADANA_LIBRARY_DIR}/lib${_lib}.so lib${_lib}.so
)
endforeach()

add_custom_command(TARGET link_prana PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink $ENV{PRAD_PATH}/database database
)
add_custom_command(TARGET link_prana PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink $ENV{PRAD_PATH}/config config
)


