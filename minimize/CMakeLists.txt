cmake_minimum_required(VERSION 3.0)

project(minimize)

set(CMAKE_CXX_STANDARD 20)

add_executable(minimize minimize.cc)
target_link_libraries(minimize PRIVATE libuni)

if(MSVC)
	target_compile_options(minimize PRIVATE /W4 /WX)
else()
	target_compile_options(minimize PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
