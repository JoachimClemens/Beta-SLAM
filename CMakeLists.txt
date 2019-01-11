CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT( bslam )

set(BSLAM_MAJOR_VERSION 1)
set(BSLAM_MINOR_VERSION 0)
set(BSLAM_PATCH_VERSION 0)
set(BSLAM_VERSION ${BSLAM_MAJOR_VERSION}.${BSLAM_MINOR_VERSION}.${BSLAM_PATCH_VERSION})
set(BSLAM_SOVERSION ${BSLAM_MAJOR_VERSION}.${BSLAM_MINOR_VERSION})

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")

# compiler settings and flags
include(CompilerSettings)


#-------------
# OpenMP
#-------------
set(BSLAM_OMP TRUE CACHE BOOL "Enable/disable OpenMP parallelization")

if(DEFINED ENV{BSLAM_OMP})
    set(BSLAM_OMP $ENV{BSLAM_OMP})
endif(DEFINED ENV{BSLAM_OMP})

if(BSLAM_OMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(BSLAM_OMP)


#-------------
# std::thread
#-------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")


#-------------
# Precompiled headers
#-------------

set(BSLAM_PRECOMPILED_HEADERS TRUE CACHE BOOL "Enable/disable precompiled headers")

if(DEFINED ENV{BSLAM_PRECOMPILED_HEADERS})
    set(BSLAM_PRECOMPILED_HEADERS $ENV{BSLAM_PRECOMPILED_HEADERS})
endif(DEFINED ENV{BSLAM_PRECOMPILED_HEADERS})

if(BSLAM_PRECOMPILED_HEADERS)
	add_definitions( -DBSLAM_PRECOMPILED_HEADERS )
    message(STATUS "Using precompiled headers.")
else()
    message(STATUS "Not using precompiled headers.")    
endif()


#-------------
# Eigen
#-------------
find_package(Eigen3 ${EIGEN3_REQUIRED_VERSION})
if(NOT EIGEN3_FOUND)
    FIND_PACKAGE(PkgConfig)
    pkg_check_modules(EIGEN3 QUIET eigen3>=${EIGEN3_REQUIRED_VERSION})
endif()
if(NOT EIGEN3_FOUND)
    message(FATAL_ERROR "Please install Eigen3.")
else()
    add_definitions( -DHAVE_EIGEN )
endif()

#-------------
# Qt5
#-------------
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Gui REQUIRED)
#include( ${QT_USE_FILE} )

#-------------
# Directories
#-------------
include_directories(
    ${INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Qt5Widgets_INCLUDE_DIRS}
    ${Qt5Gui_INCLUDE_DIRS}
)


link_directories(
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    ${PCL_LIBRARY_DIRS}
)

#--------------
# Installation
#--------------
set( BASE_DIR ${CMAKE_SOURCE_DIR} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BASE_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BASE_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BASE_DIR}/bin )

set(INSTALL_TARGETS_DEFAULT_ARGS
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)

file(MAKE_DIRECTORY ${BASE_DIR}/results)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/bslam/
	DESTINATION include/bslam/
	FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/bslam/
	DESTINATION include/bslam/
	FILES_MATCHING PATTERN "*.hpp")

add_subdirectory( src )

# Install catkin package.xml
install(FILES package.xml DESTINATION share/bslam)
  

# Finished
message("\n")
message(STATUS "Compile bslam using: make")
message(STATUS "Install bslam using: make install")
message(STATUS " (be sure to set the correct CMAKE_INSTALL_PREFIX before)")
message("\n")
