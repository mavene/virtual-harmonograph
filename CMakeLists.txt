cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

project(new_harmonograph)

list(APPEND CMAKE_CXX_FLAGS "-std=c++11")

find_package(GLFW3 REQUIRED)
find_package(GLM REQUIRED)
include_directories(${GLFW3_INCLUDE_DIR})
SET(LIBS ${GLFW3_LIBRARY})
include_directories(${GLM_INCLUDE_DIR})

#source files
file( GLOB SRCFILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)

add_subdirectory(imgui)

add_library(GLAD "src/glad.c")
set(LIBS ${LIBS} GLAD GDI32)
include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(harmonograph src/main.cpp)
target_link_libraries(harmonograph ${LIBS})
target_link_libraries(harmonograph ${GLFW3_LIBRARY})
target_link_libraries(harmonograph imgui)