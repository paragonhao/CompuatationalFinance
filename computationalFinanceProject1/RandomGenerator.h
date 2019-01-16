//
// Created by paragonhao on 15/1/19.
//
#pragma once
#ifndef COMPUTATIONALFINANCE_RANDOM_GENERATOR_H
#define COMPUTATIONALFINANCE_RANDOM_GENERATOR_H
#endif //COMPUTATIONALFINANCE_RANDOM_GENERATOR_H

class RandomGenerator {
private:
    static const int seed = 1234567890;

public:
    RandomGenerator();

    static int LGMGenerator(unsigned int m, unsigned int num);

    static long double* runif(int size);

    static long double* rbinom(int size, int n, double p);

    static long double* rexp(int size, long double *arr);

    static long double* boxmuller(long double *arr, int size);

    static long double* polarmarsaglia(long double *arr, int size);
};



