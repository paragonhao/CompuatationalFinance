//
// Created by paragonhao on 15/1/19.
//

#include <cmath>
#include "RandomGenerator.h"

using namespace std;


int RandomGenerator::LGMGenerator(unsigned int m, unsigned int num) {
    unsigned int a = pow(7, 5);
    unsigned int b = 0;
    return (a * num + b) % m ;
}

long double* RandomGenerator::runif(int size) {
    long double * arr = new long double[size];
    long double m = pow(2, 31) -1;
    arr[0] = seed;

    for(int i=0; i<=size;i++){
        arr[i+1] = RandomGenerator::LGMGenerator(m, arr[i]);
        arr[i] /= m;
    }
    return arr;
}


long double* RandomGenerator::rbinom(int size, int n, double p) {
    int uniSize = size * n;

    long double * uniformArr = RandomGenerator::runif(uniSize);
    long double * binoArr = new long double[size];

    int counter = 0;
    for(int i=0; i<size; i++){
        int sum = 0;
        // Generate bernoulli sequence and sum them up.
        for(int j=0; j<n; j++){
            sum += (uniformArr[counter] <= p)? 1 : 0;
            counter++;
        }
        binoArr[i] = sum;
    }
    return binoArr;
}

long double* RandomGenerator::rexp(int size, long double *arr){
    double lambda = 1.5;
    long double *expArr = new long double[size];

    for(int i=0; i<size; i++){
        expArr[i] = -lambda * log(arr[i]);
    }

    return expArr;
}