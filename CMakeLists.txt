cmake_minimum_required(VERSION 3.10)
project(MESHER VERSION 1.0 LANGUAGES C)
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(SOURCES
    src/main.c
    src/point.c
    src/node.c
    src/element.c
    src/export.c
)
add_executable(mesher ${SOURCES})
target_link_libraries(mesher m)
