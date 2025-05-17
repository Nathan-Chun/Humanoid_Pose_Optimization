cmake_minimum_required(VERSION 3.14)
project(Pose_Optimization_CPP)

set(CMAKE_CXX_STANDARD 17)
include_directories(${CMAKE_BINARY_DIR})
set(CMAKE_CXX_FLAGS "-O3 -ggdb -march=native -std=c++17")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


set(CASADI_InCLUDE_DIRS "/usr/local/include/casadi")
set(CASADI_LIBRARY /usr/local/lib/libcasadi.so)

find_package(CasADi REQUIRED)

include_directories(/usr/include/eigen3 /usr/local/include/casadi)

add_executable(Pose_Optimization_CPP
    src/main.cpp
    src/optimizer.cpp
    src/kinematics.cpp

)

# target_include_directories(Pose_Optimization_CPP PRIVATE usr/local/lib/libcasadi.so)
target_include_directories(Pose_Optimization_CPP PRIVATE /usr/include/eigen3 /usr/local/include/casadi ${CMAKE_SOURCE_DIR})
target_include_directories(Pose_Optimization_CPP PRIVATE ${CMAKE_SOURCE_DIR}/include/pose_optimization_cpp)
# target_link_libraries(Pose_Optimization_CPP PRIVATE${CasADi_LIBRARIES})
target_link_libraries(Pose_Optimization_CPP PRIVATE /usr/local/lib/libcasadi.so)
