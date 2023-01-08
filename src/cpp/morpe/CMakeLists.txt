cmake_minimum_required(VERSION 3.24)
project(morpe CXX)

# Pull stuff from conan
include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)

# This allows boost::stacktrace to show filenames with line numbers.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D BOOST_STACKTRACE_LINK")

# Basic setup of conan.
conan_basic_setup()

# Get all the *.cpp files
file(GLOB_RECURSE sourceFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*.cpp")

add_library(morpe STATIC ${sourceFiles} src/morpe/polynomial.cpp include/morpe.h)

# The -fPIC compiler option (gcc), requried for building libraries on linux and unix.
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

# Set the include directories
target_include_directories(${PROJECT_NAME} PUBLIC include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CONAN_INCLUDE_DIRS})

# Link to various libraries.
target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS})

# Export the project properties (including transitive dependencies) to other projects
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include
        CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
