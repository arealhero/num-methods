cmake_minimum_required(VERSION 3.0)

project(minimize)

set(CMAKE_CXX_STANDARD 20)

add_executable(minimize src/main.cc)
target_include_directories(minimize PUBLIC include)
target_compile_options(minimize PRIVATE -Wall -Wextra -Wpedantic -Werror)
