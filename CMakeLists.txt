# CMake minimum version
cmake_minimum_required (VERSION 3.20)

# Project name
project(amrFEM
        DESCRIPTION "A Hierachical structured grid FEM framework built on AMReX"
        LANGUAGES C CXX)

add_executable(amrFEM main.cpp)

#Setting the AMReX path
set(AMReX_DIR "/home/sohail/MFEM_STUFF/amrex/installdir/lib/cmake/AMReX")

# Import AMReX settings and targets
find_package(AMReX REQUIRED)
find_package(Eigen3 REQUIRED)

# Generate obj files
target_link_libraries(amrFEM PRIVATE AMReX::amrex Eigen3::Eigen)