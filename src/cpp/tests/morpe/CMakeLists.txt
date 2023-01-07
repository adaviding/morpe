cmake_minimum_required(VERSION 3.24)
project(tests_morpe)

# Pull stuff from conan
include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# Get all the *.cpp files
file(GLOB_RECURSE sourceFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*.cpp")

# Define the output type of executable.
add_executable(${PROJECT_NAME} ${sourceFiles})

# Set the include directories
target_include_directories(${PROJECT_NAME} PUBLIC include)
target_include_directories(${PROJECT_NAME} PUBLIC ${CONAN_INCLUDE_DIRS})

# Link to various libraries.
target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS})

target_include_directories(${PROJECT_NAME} PRIVATE ${morpe_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} morpe)
