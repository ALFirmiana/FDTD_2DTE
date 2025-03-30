cmake_minimum_required(VERSION 3.10)
project(FDTD_2DTE LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find HDF5 with C++ support
find_package(HDF5 REQUIRED COMPONENTS CXX)
if (NOT HDF5_FOUND)
    message(FATAL_ERROR "HDF5 library not found!")
endif()

if (NOT HDF5_CXX_FOUND)
    message(FATAL_ERROR "HDF5 C++ bindings not found! Install libhdf5-cpp-dev or similar.")
endif()

message(STATUS "HDF5 include directories: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

# Source files
set(SOURCES
    FDTD.cpp
    field.cpp
    solver.cpp
    laser.cpp
)

# Create executable
add_executable(FDTD_2DTE ${SOURCES})

# Include directories and link libraries
target_include_directories(FDTD_2DTE PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(FDTD_2DTE PRIVATE ${HDF5_LIBRARIES})

# Install target
install(TARGETS FDTD_2DTE
    RUNTIME DESTINATION bin
)
