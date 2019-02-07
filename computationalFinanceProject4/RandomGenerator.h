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

    static long double* rexp(long double *arr, int size, double lambda);

    static long double* boxmuller(long double *arr, int size);

    static long double * boxmullerWithHaltonSeq(long double *base1, long double * base2,int size);

    static long double* polarmarsaglia(long double *arr, int size);

    static long double*  bivariateNormalX(long double *z1, int size);

    static long double*  bivariateNormalY(long double *z1, long double *z2, double rho, int size);

    static long double* wienerProcess(double t, int size, int seed);

    static long double * getHaltonSequence(int base, int size);
};



