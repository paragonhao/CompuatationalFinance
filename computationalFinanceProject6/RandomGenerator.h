//
// Created by paragonhao on 15/1/19.
//
#pragma once
#ifndef COMPUTATIONALFINANCE_RANDOM_GENERATOR_H
#define COMPUTATIONALFINANCE_RANDOM_GENERATOR_H


class RandomGenerator {

public:

    static const int seed = 1234567890;

    RandomGenerator();

    static int LGMGenerator(unsigned int m, unsigned int num);

    static double* runif(int size, int seed);

    static double* rbinom(int size, int n, double p, int seed);

    static double* rexp(double *arr, int size, double lambda);

    static double* boxmuller(double *arr, int size);

    static double * boxmullerWithHaltonSeq(double *base1, double * base2,int size);

    static double* polarmarsaglia(double *arr, int size);

    static double*  bivariateNormalX(double *z1, int size);

    static double*  bivariateNormalY(double *z1, double *z2, double rho, int size);

    static double* wienerProcess(double t, int size, int seed);

    static double * getHaltonSequence(int base, int size);
};
#endif //COMPUTATIONALFINANCE_RANDOM_GENERATOR_H


