cmake_minimum_required(VERSION 3.29)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(deps/googletest-main)
project(Chemilator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)
include_directories(deps/googletest-main/googlemock/include)
include_directories(deps/googletest-main/googletest/include)

FILE(GLOB SRC_CONTENTS "src/*.cpp")

add_library(chemilator STATIC ${SRC_CONTENTS} include/chm/sqlite3.h src/sqlite3.c)
target_link_libraries(chemilator)

add_executable(main_test test/main_test.cpp)
target_link_libraries(main_test chemilator gtest gtest_main)