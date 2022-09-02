cmake_minimum_required(VERSION 3.10)
project(MyOdeSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED NO_MODULE)


set(SOURCES
  src/main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME}
  PRIVATE "${PROJECT_SOURCE_DIR}/include"
  PUBLIC ${EIGEN3_INCLUDE_DIR}
)

