# from within embuild directory
# emcmake cmake .. && emmake make

cmake_minimum_required(VERSION 3.29)
project(Chemilator)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Ensure the output directory exists
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/embuild)

include_directories(include)

FILE(GLOB SRC_CONTENTS "src/*.cpp")

add_library(chemilator STATIC ${SRC_CONTENTS})

add_executable(eqnBalance impl/eqnBalance.cpp)
target_link_libraries(eqnBalance chemilator)


# Check if the compiler is Emscripten
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Emscripten")
    message("EMSCRIPTEN compiler is in use")

    set_target_properties(eqnBalance PROPERTIES
        LINK_FLAGS "-s WASM=1 -s MODULARIZE=1 -s EXPORT_ALL=1 -o ${CMAKE_SOURCE_DIR}/pages/wasm/eqnBalance.js --shell-file ${CMAKE_SOURCE_DIR}/shell.html"
    )
else()
    # message(FATAL_ERROR "This project requires Emscripten as the compiler")
endif()
