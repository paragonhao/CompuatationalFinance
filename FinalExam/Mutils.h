//
// Created by paragonhao on 16/1/19.
//

#ifndef COMPUTATIONALFINANCE_MUTILS_H
#define COMPUTATIONALFINANCE_MUTILS_H

#include <string>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class Mutils {
public:
    static double Mean(double arr[], int n);
    static double StDev(double arr[], int n);
    static double Corr(double *x, double *y, int size);
    static void WriteArrayToCSV( double *arr, int n, const string filename);
    static void WriteToCSV2DMatrix( double arr[][31], int col,int row, const string filename);
    static double pnorm(double x);
    static double* runif(int size, int seed);
    static double Cov(double *x, double *y, int size);
    static double * MatrixMultiply(double *arr, int size, double t);
    static double * MatrixAddition(double *arr, int size, double num);
    static double max(const double & num1, const double & num2);
    static double min(const double & num1, const double & num2);
    static double arrayElementWiseMultiply(double * v1, double * v2, int size);
    static void cumSum(MatrixXd &matrix, MatrixXd &cumSumMatrix, int nSteps, int nSim);
    static void generateWProcessMat(int nSim, int nSteps, double delta_t, MatrixXd &mat, int seed);
    static void generatePricePath(int nSteps, int nSim, double T, double s0, double r, double sigma, MatrixXd &priceProcess);

};


#endif //COMPUTATIONALFINANCE_MUTILS_H
