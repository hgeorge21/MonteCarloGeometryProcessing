cmake_minimum_required(VERSION 3.1)
project(MonteCarlo-PDE)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) 

if(APPLE)
    set(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
    set(OPENMP_INCLUDES "/usr/local/opt/llvm/include")
    set(OpenMP_C "${CMAKE_C_COMPILER}")
    set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
    set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})

    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3")
endif()
if(UNIX AND NOT APPLE)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
endif()

### libIGL options: choose between header only and compiled static library
set(LIBIGL ${PROJECT_SOURCE_DIR}/libigl)
set(LIBIGL_INCLUDE_DIR ${LIBIGL}/include)
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

# Can install a fuller version of eigen @ https://gitlab.com/libeigen/eigen.git
set(EIGENROOT "${PROJECT_SOURCE_DIR}/libigl/external/eigen")
include_directories("${EIGENROOT}")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
if(USE_SOLUTION)
  file(GLOB SRCFILES solution/*.cpp)
else()
  file(GLOB SRCFILES src/*.cpp)
endif()

add_executable(${PROJECT_NAME} "main.cpp" ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw )
