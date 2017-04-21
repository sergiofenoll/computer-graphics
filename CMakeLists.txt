cmake_minimum_required(VERSION 3.6)
# Project name
project(CompGraphX)
# For vim autocomplete
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# -std=c++11
set(CMAKE_CXX_STANDARD 11)
# Source files
set(SOURCE_FILES
        src/easy_image.cc
        src/easy_image.hh
        src/engine.cc
        src/ini_configuration.cc
        src/ini_configuration.hh
        src/l_parser.cc
        src/l_parser.hh
        src/Line2D.cc
        src/Line2D.hh
        src/3DFigures.cc
        src/3DFigures.hh
        src/vector.cc
        src/vector.hh)
# Single executable
add_executable(engine ${SOURCE_FILES})
# Create tgz archive
add_custom_target(tgz COMMAND
    ${CMAKE_COMMAND} -E tar "cfvz" "../engine.tgz"
       "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
       "${CMAKE_CURRENT_SOURCE_DIR}/cmake-build-debug/Makefile"
       "${CMAKE_CURRENT_SOURCE_DIR}/src/")
