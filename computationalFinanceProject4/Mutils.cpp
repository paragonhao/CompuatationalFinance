//
// Created by paragonhao on 16/1/19.
//

#include "Mutils.h"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double Mutils::Mean(long double arr[], int n)
{
    double mean = 0.0;
    for (int idx = 0; idx < n; idx++)
    {
        mean += arr[idx];
    }
    mean /= static_cast<double>(n);
    return mean;
}


double Mutils::StDev(long double arr[], int n)
{
    double mean = Mutils::Mean(arr, n);
    long double variance = 0.0;
    for (size_t idx = 0; idx < n; idx++)
    {
        double temp = arr[idx] - mean;
        variance += temp*temp;
    }

    // Compute sample variance using Bessel's correction (see http://en.wikipedia.org/wiki/Bessel%27s_correction)
    variance /= static_cast<long double>(n) - (n == 1 ? 0.0 : 1.0);

    // Standard deviation is square root of variance
    return std::sqrt(variance);
}

double Mutils::Corr(long double *x, long double *y, int size)
{
    double meanX = Mutils::Mean(x, size);
    double meanY = Mutils::Mean(y, size);
    double numerator = 0.0;
    double denominatorX = 0.0;
    double denominatorY = 0.0;

    for(int i = 0; i< size; i++){
        numerator += (x[i] - meanX) * (y[i] - meanY);
        denominatorX += (x[i] - meanX) * (x[i] - meanX);
        denominatorY += (y[i] - meanY) * (y[i] - meanY);
    }

    numerator /= (size-1);
    denominatorX = sqrt(denominatorX /(size - 1));
    denominatorY = sqrt(denominatorY /(size -1));

    return numerator/(denominatorY* denominatorX);
}

double Mutils::Cov(long double *x, long double *y, int size){
    double meanX = Mutils::Mean(x, size);
    double meanY = Mutils::Mean(y, size);
    double sse = 0 ;
    for(int i = 0; i< size; i++){
        sse += (x[i] - meanX) * (y[i] - meanY);
    }
    return sse/(size-1);
}



void Mutils::WriteArrayToCSV( double *arr, int n, const string filename){
    ofstream file;
    file.open(filename);
    for(int i=0; i < n; i++){
        file << arr[i] << ",\n";
    }
    file.close();
}


void Mutils::WriteToCSV2DMatrix(double arr[][31], int row,int col, const string filename){
    ofstream file;
    file.open(filename);

    for(int i=0; i<row;i++){
        for(int j=0; j<col; j++){
            file << arr[j][i] << ",";
        }
        file << "\n";
    }

    file.close();
}

long double * Mutils::MatrixMultiply(long double *arr, int size, double t){

    for(int i=0; i< size; i++){
        arr[i] *= t;
    }

    return arr;
}

long double * Mutils::MatrixAddition(long double *arr, int size, long double num){

    for(int i=0; i< size; i++){
        arr[i] += num;
    }

    return arr;
}

long double Mutils::pnorm(double x)
{
    const long double d1 = 0.0498673470;
    const long double d2 = 0.0211410061;
    const long double d3 = 0.0032776263;
    const long double d4 = 0.0000380036;
    const long double d5 = 0.0000488906;
    const long double d6 = 0.0000053830;

    long double probability = 0.0;
    long double absX = (x >=0)?x:-x;
    long double temp =  1 + (d1 * absX) + (d2 * absX * absX) + (d3 * absX * absX * absX) + (d4 * pow(absX,4))
            + (d5 * pow(absX,5)) + (d6 * pow(absX, 6));
    probability = 1 - 0.5 * pow(temp, -16);

    return  (x >= 0)? probability : (1- probability);
}


/* Vector element wise multiplication
 * s
 *
 * @param:
 * double *: first array
 * double *: seconed array
 * int : size of the array to check
 * */
long double Mutils::arrayElementWiseMultiply(long double * v1, long double * v2, int NumBits){
    long double sum = 0;
    for(int i = 0; i<NumBits; i++){
        sum += v1[i] * v2[i];
    }
    return sum;
}
