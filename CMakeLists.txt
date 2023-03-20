cmake_minimum_required(VERSION 3.26)

project(
  ModernCMakeExample
  VERSION 0.1
  DESCRIPTION "An example project with CMake"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(SOURCES 
    ./src/main.c++
    ./src/Field.c++
    ./src/Mesh.c++
    ./src/Model.c++)

set(HEADERS
    ./include/Field.h
    ./include/Mesh.h
    ./include/Model.h)

add_executable(main ${SOURCES} ${HEADERS})

target_include_directories(
    main
    PRIVATE
    ./include/
)

