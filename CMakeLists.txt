cmake_minimum_required (VERSION 3.20)

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

project (num-methods)

set (CMAKE_CXX_STANDARD 20 REQUIRED)
set (CMAKE_CXX_FLAGS_DEBUG "-g")
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

if (USE_CLANG_TIDY)
  set (CMAKE_CXX_CLANG_TIDY clang-tidy ;--format-style=file;--use-color)
endif ()

add_subdirectory (vendors)
include_directories (SYSTEM vendors/eigen)

include (compile_options)

include_directories (lib)
add_subdirectory (lib)

add_subdirectory (04)
add_subdirectory (05)
add_subdirectory (tests)

