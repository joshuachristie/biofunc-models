cmake_minimum_required(VERSION 3.20.2)

project(QEF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -O3")
file(GLOB SOURCES
    *.cpp
    *.h
    include/*.cc
    include/*.h
)

add_executable(QEF ${SOURCES})

find_package(Protobuf REQUIRED)
target_link_libraries(QEF ${Protobuf_LIBRARIES})

