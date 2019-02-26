cmake_minimum_required(VERSION 3.10)
project(computationalFinanceProject7)

set(CMAKE_CXX_STANDARD 17)
set(EIGEN_DIR "/Users/paragonhao/Documents/dev/eigen")
include_directories(${EIGEN_DIR})

add_executable(computationalFinanceProject7 main.cpp RandomGenerator.h RandomGenerator.cpp Mutils.cpp Mutils.h OptionPricing.cpp OptionPricing.h JDDefaultOption.cpp JDDefaultOption.h DifferenceMethod.cpp DifferenceMethod.h)