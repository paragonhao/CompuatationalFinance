//
// Created by paragonhao on 15/1/19.
//

#include <cmath>
#include <iostream>
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

// Simulate Standard Normal Distribution Using Box-Muller
long double* RandomGenerator::boxmuller(long double *arr, int size){

    long double* normArr = new long double[size];

    for(int i=0; i< size; i+=2){
        // calculate z1
        normArr[i] = sqrt(-2 * log(arr[i])) * cos(2 * M_PI * arr[i+1]);

        // calculate z2
        normArr[i+1] = sqrt(-2 * log(arr[i])) * sin(2 * M_PI * arr[i+1]);
    }

    return normArr;
}

// Simulate Standard Normal Distribution Using Polar-Marsaglia
// the first item in the array contains the number of random variables available
long double* RandomGenerator::polarmarsaglia(long double *arr, int size){

    long double* normArr = new long double[size];
    int arrSize = 1;
    long double v1 = 0;
    long double v2 = 0;
    long double w = 0;


    for(int i=0; i<size; i+=2){
        v1 = 2 * arr[i] - 1;
        v2 = 2 * arr[i+1] -1;
        w = v1 * v1  + v2 * v2; // using pow() would significantly increase the execution time.

        if(w <= 1.0 ){
            normArr[arrSize] = v1 * sqrt(-2 * log(w) / w);
            normArr[arrSize+1] = v2 * sqrt(-2 * log(w) / w);
            arrSize +=2;
        }
    }

    //Put the total number of random variables generated in the first cell.
    normArr[0] = arrSize -1;

    return normArr;
}