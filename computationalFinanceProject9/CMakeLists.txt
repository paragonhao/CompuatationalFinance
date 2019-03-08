cmake_minimum_required(VERSION 3.10)
project(computationalFinanceProject9)

set(CMAKE_CXX_STANDARD 17)

set(EIGEN_DIR "/Users/paragonhao/Documents/dev/eigen")
include_directories(${EIGEN_DIR})

add_executable(computationalFinanceProject9 main.cpp MortgageBackedSecurities.cpp MortgageBackedSecurities.h
        FixedIncome.cpp FixedIncome.h Mutils.h Mutils.cpp)