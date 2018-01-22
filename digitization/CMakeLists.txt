#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PRadDig)

#----------------------------------------------------------------------------
# Tell CMake where to find the ROOT installation
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED)

#---------------------------------------------------------------------------
# Locate EVIO library and headers
#
if(DEFINED ENV{ET_INC} AND DEFINED ENV{ET_LIB})
  set(EVIO_INCLUDE_DIR $ENV{ET_INC} CACHE PATH "Path to evio headers")
  set(EVIO_LIBRARY_DIR $ENV{ET_LIB} CACHE PATH "Path to evio libraries")
elseif(DEFINED ENV{EVIODIR})
  set(EVIO_INCLUDE_DIR $ENV{EVIODIR}/include CACHE PATH "Path to evio headers")
  set(EVIO_LIBRARY_DIR $ENV{EVIODIR}/lib CACHE PATH "Path to evio libraries")
else()
  set(EVIO_INCLUDE_DIR EVIO_INCLUDE_DIR_NOT_FOUND CACHE PATH "Path to evio headers")
  set(EVIO_LIBRARY_DIR EVIO_INCLUDE_DIR_NOT_FOUND CACHE PATH "Path to evio libraries")
endif()
set(EVIO_LIBRARIES evio)

#---------------------------------------------------------------------------
# Locate pradanalyzer library and headers
#
if(DEFINED ENV{PRAD_PATH})
  set(PRADANA_INCLUDE_DIR $ENV{PRAD_PATH}/include CACHE PATH "Path to prad analyzer headers")
  set(PRADANA_LIBRARY_DIR $ENV{PRAD_PATH}/lib CACHE PATH "Path to prad analyzer libraries")
elseif(DEFINED ENV{PRADANADIR})
  set(PRADANA_INCLUDE_DIR $ENV{PRADANADIR}/include CACHE PATH "Path to prad analyzer headers")
  set(PRADANA_LIBRARY_DIR $ENV{PRADANADIR}/lib CACHE PATH "Path to prad analyzer libraries")
else()
  set(PRADANA_INCLUDE_DIR PRADANA_INCLUDE_DIR_NOT_FOUND CACHE PATH "Path to prad analyzer headers")
  set(PRADANA_LIBRARY_DIR PRADANA_INCLUDE_DIR_NOT_FOUND CACHE PATH "Path to prad analyzer libraries")
endif()
set(PRADANA_LIBRARIES prana prconf)

#----------------------------------------------------------------------------
# Setup include and lib directory for this project
#
include_directories(${ROOT_INCLUDE_DIRS})
if(NOT ${EVIO_INCLUDE_DIR} STREQUAL "EVIO_INCLUDE_DIR_NOT_FOUND" OR ${EVIO_LIBRARY_DIR} STREQUAL "EVIO_LIBRARY_DIR_NOT_FOUND")
  include_directories(${EVIO_INCLUDE_DIR})
  link_directories(${EVIO_LIBRARY_DIR})
else()
  message(FATAL_ERROR "Evio headers and libraries not found")
endif()
if(NOT ${PRADANA_INCLUDE_DIR} STREQUAL "PRADANA_INCLUDE_DIR_NOT_FOUND" OR ${PRADANA_LIBRARY_DIR} STREQUAL "PRADANA_LIBRARY_DIR_NOT_FOUND")
  include_directories(${PRADANA_INCLUDE_DIR})
  link_directories(${PRADANA_LIBRARY_DIR})
else()
  message(FATAL_ERROR "PRad Analyzer headers and libraries not found")
endif()
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB PRADDIG_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB PRADDIG_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries
#
add_executable(PRadDig PRadDig.cc ${PRADDIG_SOURCES} ${PRADDIG_HEADERS})
target_link_libraries(PRadDig ${ROOT_LIBRARIES} ${EVIO_LIBRARIES} ${PRADANA_LIBRARIES})
add_executable(PRadRec PRadRec.cc)
target_link_libraries(PRadRec ${ROOT_LIBRARIES} ${PRADANA_LIBRARIES})

#----------------------------------------------------------------------------
# Copy the executable to the project directory
#
add_custom_target(copy_exes ALL
  DEPENDS PRadDig PRadRec
  COMMENT "Copy executable"
)
add_custom_command(TARGET copy_exes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy PRadDig ..
)
add_custom_command(TARGET copy_exes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy PRadRec ..
)
add_custom_command(TARGET copy_exes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${PRADANA_LIBRARY_DIR}/libprana.so ../libprana.so
)
