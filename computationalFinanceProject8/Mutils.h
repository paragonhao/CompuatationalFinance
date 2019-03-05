//
// Created by paragonhao on 16/1/19.
//

#ifndef COMPUTATIONALFINANCE_MUTILS_H
#define COMPUTATIONALFINANCE_MUTILS_H

#include <string>

using namespace std;

class Mutils {
public:
    static double Mean(double arr[], int n);
    static double StDev(double arr[], int n);
    static double Corr(double *x, double *y, int size);
    static void WriteArrayToCSV( double *arr, int n, const string filename);
    static void WriteToCSV2DMatrix( double arr[][31], int col,int row, const string filename);
    static double pnorm(double x);
    static double Cov(double *x, double *y, int size);
    static double * MatrixMultiply(double *arr, int size, double t);
    static double * MatrixAddition(double *arr, int size, double num);
    static double max(const double & num1, const double & num2);
    static double arrayElementWiseMultiply(double * v1, double * v2, int size);
};


#endif //COMPUTATIONALFINANCE_MUTILS_H
