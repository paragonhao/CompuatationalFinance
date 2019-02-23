cmake_minimum_required(VERSION 3.10)
project(ComputationalFinanceProject6)

set(CMAKE_CXX_STANDARD 17)

add_executable(ComputationalFinanceProject6 main.cpp RandomGenerator.h RandomGenerator.cpp Mutils.cpp Mutils.h OptionPricing.cpp OptionPricing.h JDDefaultOption.cpp JDDefaultOption.h)