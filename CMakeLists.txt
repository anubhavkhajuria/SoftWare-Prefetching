cmake_minimum_required(VERSION 3.12)
project(PrefetchIMA_Pass)


set(CMAKE_CXX_STANDARD 17)


find_package(LLVM REQUIRED CONFIG)
include(AddLLVM)

add_definitions(${LLVM_DEFINITIONS})
include_directories(${LLVM_INCLUDE_DIRS})
link_directories(${LLVM_LIBRARY_DIRS})

if(APPLE)
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif()

add_subdirectory(pass)

