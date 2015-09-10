find_package(PkgConfig REQUIRED)

cmake_minimum_required(VERSION 3.3)
project(nsgt)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-D_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=undefined")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -g -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls -fsanitize=undefined")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

pkg_check_modules(LIBSNDFILE REQUIRED sndfile) #package libsndfile-dev
include_directories(${LIBSNDFILE_INCLUDE_DIRS})
link_directories(${LIBSNDFILE_LIBRARY_DIRS})

pkg_check_modules(FFTW REQUIRED fftw3) #package libfftw3-dev
include_directories(${FFTW_INCLUDE_DIRS})
link_directories(${FFTW_LIBRARY_DIRS})

set(SOURCE_FILES main.cpp nsgt.c nsgt.h examples/example_spectrogram.h examples/example_nsgt_error.h examples/example_utils.h examples/example_nsgt_stream_error.h)
add_executable(nsgt ${SOURCE_FILES})

target_link_libraries(nsgt ${LIBSNDFILE_LIBRARIES} ${FFTW_LIBRARIES})