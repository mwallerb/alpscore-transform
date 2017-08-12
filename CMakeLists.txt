cmake_minimum_required(VERSION 2.8.10)

project (CTINT)

# see https://blog.kitware.com/upcoming-in-cmake-2-8-12-osx-rpath-support/
set(CMAKE_MACOSX_RPATH 1)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

find_package(Eigen3 REQUIRED)
find_package(FFTW3)
find_package(GTest)
find_package(SWIG)
find_package(PythonLibs)
find_package(NumPy)

if(SWIG_FOUND AND PYTHONLIBS_FOUND AND PYTHON_NUMPY_FOUND)
    set(SWIG_NUMPY_FOUND TRUE)
else()
    set(SWIG_NUMPY_FOUND FALSE)
endif()

option(WITH_FFTW3 "Build with FFTW3 support" ${FFTW3_FOUND})
option(WITH_SWIG "Create SWIG wrappers" ${SWIG_NUMPY_FOUND})
option(WITH_GTEST "Build gtest tests" ${GTEST_FOUND})

include_directories(${FFTW3_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

if(WITH_FFTW3)
    add_definitions(-DHAVE_FFTW=1)
else()
    add_definitions(-DHAVE_FFTW=0)
endif()

if(NOT CMAKE_BUILD_TYPE)
    message(WARNING "Preparing DEFAULT build (you want Release or Debug)...")
    add_definitions(-DMODE_DEBUG=0)
    add_definitions(-DALPS_GF_DEBUG)
elseif(CMAKE_BUILD_TYPE MATCHES ^Release|RelWithDebInfo$)
    message("Preparing RELEASE build ...")
    add_definitions(-DMODE_DEBUG=0)
    add_definitions(-DBOOST_DISABLE_ASSERTS)
    add_definitions(-DNDEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -flto")
elseif(CMAKE_BUILD_TYPE STREQUAL Debug)
    message("Preparing DEBUG build ...")
    add_definitions(-DMODE_DEBUG=1)
    add_definitions(-DALPS_GF_DEBUG)
else()
    message(FATAL_ERROR "Unknown build type: ${CMAKE_BUILD_TYPE}" )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

include_directories(include)

add_library(transform-lib SHARED
            src/fftw.cpp src/fourier.cpp src/nonuniform.cpp)
set_target_properties(transform-lib PROPERTIES
                      OUTPUT_NAME alps-transform
                      INTERPROCEDURAL_OPTIMIZATION True
                      )

target_link_libraries(transform-lib
                      ${FFTW3_LIBRARIES}
                      )

# See <http://stackoverflow.com/questions/16074598>
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH
         "Default install path" FORCE)
endif()

install(DIRECTORY include/
        DESTINATION include/
        FILES_MATCHING PATTERN "*.hpp"
        )
install(TARGETS transform-lib
        RUNTIME DESTINATION "bin/"
        LIBRARY DESTINATION "lib/"
        )

if(PYTHONLIBS_FOUND OR WITH_SWIG)
    string(REGEX MATCH "^[0-9]+\\.[0-9]+" PYTHON_VERSION
           "${PYTHONLIBS_VERSION_STRING}")
    set(PYTHON_INSTALL_DIR "lib/python${PYTHON_VERSION}/site-packages"
        CACHE PATH "Path for installing the python packages")

    install(DIRECTORY python/
            DESTINATION ${PYTHON_INSTALL_DIR}/alps
            FILES_MATCHING PATTERN "*.py"
            )
endif()

if(WITH_SWIG)
    include(${SWIG_USE_FILE})

    set(CMAKE_SWIG_FLAGS "")

    set_source_files_properties(swig/transform.i PROPERTIES
                                CPLUSPLUS ON
                                #SWIG_FLAGS "-includeall"
                                )
    #file(GLOB SWIG_MODULE_ctint_EXTRA_DEPS swig/*.inc.i)

    swig_add_module(transform python swig/transform.i)
    swig_link_libraries(transform
                        transform-lib
                        ${PYTHON_LIBRARIES}
                        ${FFTW3_LIBRARIES}
                        )

    set_source_files_properties(${swig_generated_file_fullname} PROPERTIES
                                COMPILE_FLAGS "-Wno-long-long -Wno-extra")

    set_target_properties(_transform PROPERTIES
                          INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/"
                          LINK_FLAGS "-Wno-long-long -Wno-extra")
    install(TARGETS _transform
            LIBRARY DESTINATION ${PYTHON_INSTALL_DIR}/alps
            )
    install(FILES ${CMAKE_BINARY_DIR}/transform.py
            DESTINATION ${PYTHON_INSTALL_DIR}/alps
            )
endif()

if(WITH_GTEST)
    enable_testing()
    add_subdirectory(test)
endif()
