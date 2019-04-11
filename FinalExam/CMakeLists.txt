cmake_minimum_required(VERSION 3.10)
project(FinalExam)

set(CMAKE_CXX_STANDARD 17)

set(EIGEN_DIR "/Users/paragonhao/Documents/dev/eigen")
include_directories(${EIGEN_DIR})

add_executable(FinalExam main.cpp Final_1.h Final_2.h Final_3.h Final_4.h Final_5.h Mutils.h Mutils.cpp)