cmake_minimum_required(VERSION 3.10)
project(computationalfinance)

set(CMAKE_CXX_STANDARD 17)

add_executable(computationalfinance main.cpp RandomGenerator.h RandomGenerator.cpp Mutils.cpp Mutils.h OptionPricing.cpp OptionPricing.h)