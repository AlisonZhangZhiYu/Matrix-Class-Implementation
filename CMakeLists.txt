cmake_minimum_required(VERSION 3.28)
project(MatrixClass)

set(CMAKE_CXX_STANDARD 17)

add_executable(MatrixClass main.cpp
        Matrix.cpp
        Matrix.h)
