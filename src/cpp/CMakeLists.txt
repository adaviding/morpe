cmake_minimum_required(VERSION 3.24)
project(morpe_cpp CXX)

set(CMAKE_CXX_STANDARD 20)

# Set the default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)  # Possible values are:  Debug, Release, RelWithDebInfo and MinSizeRel
endif()

# Compiler flags
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # using Clang
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wno-unknown-pragmas -Wno-unused-label -Wno-attributes")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-unknown-pragmas -Wno-unused-label -Wno-attributes")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # using GCC
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wno-unknown-pragmas -Wno-unused-label -Wno-attributes")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-unknown-pragmas -Wno-unused-label -Wno-attributes")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # using Visual Studio C++
endif()

# Remove 'lib' prefix for Windows libraries
if (WIN32)
    # This affects the names of the libraries that are being pulled in (as dependencies).
    set(CMAKE_SHARED_LIBRARY_PREFIX "")
    set(CMAKE_STATIC_LIBRARY_PREFIX "")
endif ()

# Sub-projects:  The order matters.
add_subdirectory(morpe)
add_subdirectory(tests)
