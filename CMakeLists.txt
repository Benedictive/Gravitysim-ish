cmake_minimum_required(VERSION 3.12)

# set the project name and version
project(GravitySim
        VERSION 1.0
        LANGUAGES CXX CUDA)

find_package(SDL2 REQUIRED CONFIG REQUIRED COMPONENTS SDL2)
find_package(SDL2 REQUIRED CONFIG COMPONENTS SDL2main)

add_subdirectory(extern/SimpleCanvas)

add_subdirectory(src)

add_subdirectory(apps)