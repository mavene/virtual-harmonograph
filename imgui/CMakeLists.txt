
file(GLOB_RECURSE source *.cpp)
file(GLOB_RECURSE headerfile *.h)
add_library(imgui ${source})
target_include_directories(imgui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})