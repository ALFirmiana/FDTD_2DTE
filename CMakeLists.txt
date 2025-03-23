cmake_minimum_required(VERSION 3.10)
project(FDTD_2DTE LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(HDF5 REQUIRED)

# Include directories
include_directories(${HDF5_INCLUDE_DIRS})

# Source files
set(SOURCES
    FDTD.cpp
    field.cpp
    solver.cpp
)

# Create executable
add_executable(FDTD_2DTE ${SOURCES})

# Link libraries
target_link_libraries(FDTD_2DTE PRIVATE ${HDF5_LIBRARIES})

# Install target
install(TARGETS FDTD_2DTE
    RUNTIME DESTINATION bin
)

