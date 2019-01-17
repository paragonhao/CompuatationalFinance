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
    static void WriteToCSV(long double *arr, int n, const string filename);
};


#endif //COMPUTATIONALFINANCE_MUTILS_H
