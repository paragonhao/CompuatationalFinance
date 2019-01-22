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

double Mutils::Cov(long double x[], long double y[], int size)
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



void Mutils::WriteToCSV(long double *arr, int n, const string filename){
    ofstream file;
    file.open(filename);
    for(int i=0; i < n; i++){
        file << arr[i] << ",\n";
    }
    file.close();
}