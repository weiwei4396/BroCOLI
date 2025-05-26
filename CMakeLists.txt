#CMake版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(BroCOLI)

# 设置 C++ 标准为 C++11
set(CMAKE_CXX_STANDARD 11)

# 设置线程支持
find_package(Threads REQUIRED)

# 添加 Eigen 头文件路径
# 假设 Eigen 库已经安装到 /data/workdir/panw/C_code/source/ 文件夹下
include_directories(${CMAKE_SOURCE_DIR}/external/Eigen)

# 添加源文件
set(SOURCES BroCOLI_SC_v2_column_FullLength.cpp)

# 创建可执行文件
add_executable(BroCOLI ${SOURCES})

# 链接线程库
target_link_libraries(BroCOLI Threads::Threads)

