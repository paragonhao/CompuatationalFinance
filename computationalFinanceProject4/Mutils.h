//
// Created by paragonhao on 16/1/19.
//

#ifndef COMPUTATIONALFINANCE_MUTILS_H
#define COMPUTATIONALFINANCE_MUTILS_H

#include <string>

using namespace std;

class Mutils {
public:
    static double Mean(long double arr[], int n);
    static double StDev(long double arr[], int n);
    static double Corr(long double *x, long double *y, int size);
    static void WriteArrayToCSV( double *arr, int n, const string filename);
    static void WriteToCSV2DMatrix( double arr[][31], int col,int row, const string filename);
    static long double pnorm(double x);
    static double Cov(long double *x, long double *y, int size);
    static long double * MatrixMultiply(long double *arr, int size, double sqrtT);
    static long double * MatrixAddition(long double *arr, int size, long double num);
    static long double * getHaltonSequence(int base, int size);
    static long double arrayElementWiseMultiply(long double * v1, long double * v2, int size);
};


#endif //COMPUTATIONALFINANCE_MUTILS_H
