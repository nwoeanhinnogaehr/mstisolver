cmake_minimum_required(VERSION 3.1...3.27)

project(MSTISolver)
add_compile_options(
  -std=c++17
  -march=native
  -mtune=native
)
include_directories(${CMAKE_SOURCE_DIR}/3rdparty/include)
add_executable(mstisolver src/main.cpp)
