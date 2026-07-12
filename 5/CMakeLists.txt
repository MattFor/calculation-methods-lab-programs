cmake_minimum_required(VERSION 4.2)
project(Program5)

set(CMAKE_CXX_STANDARD 26)

include(FetchContent)

FetchContent_Declare(
        hprint
        GIT_REPOSITORY https://github.com/MattFor/hprint.git
        GIT_TAG v1.1.0
)

FetchContent_MakeAvailable(hprint)

add_executable(Program5 main.cpp)

target_link_libraries(Program5 PRIVATE hp::hprint)