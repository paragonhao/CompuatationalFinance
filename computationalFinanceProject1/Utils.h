//
// Created by paragonhao on 15/1/19.
//
#pragma once
#ifndef COMPUTATIONALFINANCE_STATSUTILS_H
#define COMPUTATIONALFINANCE_STATSUTILS_H

#include <cmath>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;

template <typename T>
double Mean(const T arr[], size_t n)
{
    double mean = 0.0;
    for (size_t idx = 0; idx < n; idx++)
    {
        mean += arr[idx];
    }
    mean /= static_cast<double>(n);
    return mean;
}

template <typename T>
double StDev(const T arr[], size_t n)
{
    double mean = Mean(arr, n);
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

template <typename T>
void WriteToCSV(const T *arr, size_t n, const string filename){
    ofstream file;
    file.open(filename);
    for(int i=0; i < n; i++){
        file << arr[i] << ",\n";
    }
    file.close();
}

#endif //COMPUTATIONALFINANCE_STATSUTILS_H
