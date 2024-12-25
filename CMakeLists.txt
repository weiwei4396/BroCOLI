#CMake版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(BroCOLI VERSION 1.0 LANGUAGES CXX)
# 设置 C++ 标准为 C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# 设置线程支持
find_package(Threads REQUIRED)

# 指定 src 文件夹为包含目录（用于头文件）
include_directories(${CMAKE_SOURCE_DIR}/src)

# 创建可执行文件
add_executable(BroCOLI_bulk src/BroCOLI_Bulk.cpp)
target_link_libraries(BroCOLI_bulk PRIVATE Threads::Threads)

add_executable(BroCOLI_sc src/BroCOLI_SingleCell.cpp)
target_link_libraries(BroCOLI_sc PRIVATE Threads::Threads)

