//
// Created by paragonhao on 15/1/19.
//
#pragma once
#ifndef COMPUTATIONALFINANCE_RANDOM_GENERATOR_H
#define COMPUTATIONALFINANCE_RANDOM_GENERATOR_H


class RandomGenerator {
private:
    static const int seed = 1234567890;

public:
    RandomGenerator(int = 1234);

    static int LGMGenerator(unsigned int m, unsigned int num);

    static long double* runif(int size);

    static long double* rbinom(int size, int n, double p);

    static long double* rexp(int size, long double *arr);
};


#endif //COMPUTATIONALFINANCE_RANDOM_GENERATOR_H
