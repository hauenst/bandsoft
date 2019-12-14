# CMakeLists for tools
project(tools)

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/dependencies/cmake_modules")
find_package(ROOT REQUIRED)
find_package(LZ4 REQUIRED)
add_definitions(-D__LZ4__)

# Add directories
add_subdirectory(dependencies)
add_subdirectory(tools)
add_subdirectory(calibration)
add_subdirectory(analysis)
#add_subdirectory(simulation)
