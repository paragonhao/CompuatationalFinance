//
// Created by paragonhao on 15/1/19.
//
#pragma once
#ifndef COMPUTATIONALFINANCE_RANDOM_GENERATOR_H
#define COMPUTATIONALFINANCE_RANDOM_GENERATOR_H
#endif //COMPUTATIONALFINANCE_RANDOM_GENERATOR_H

class RandomGenerator {

public:

    static const int seed = 1234567890;

    RandomGenerator();

    static int LGMGenerator(unsigned int m, unsigned int num);

    static long double* runif(int size, int seed);

    static long double* rbinom(int size, int n, double p, int seed);

    static long double* rexp(int size, long double *arr);

    static long double* boxmuller(long double *arr, int size);

    static long double* polarmarsaglia(long double *arr, int size);

    static long double*  bivariateNormalX(long double *z1, double size);

    static long double*  bivariateNormalY(long double *z1, long double *z2, double rho, double size);

    static long double* wienerProcess(double t, int size, int seed);
};



