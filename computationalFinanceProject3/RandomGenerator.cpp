//
// Created by paragonhao on 15/1/19.
//

#include <cmath>
#include <iostream>
#include "RandomGenerator.h"
#include "Mutils.h"
using namespace std;


int RandomGenerator::LGMGenerator(unsigned int m, unsigned int num) {
    unsigned int a = pow(7, 5);
    unsigned int b = 0;
    return (a * num + b) % m ;
}


/* Generate Uniform Distribution using LGM generator
 *
 * @param:
 * int size: array size to be generated
 * int seed: seed which is used to generate standard normal
 * */
long double* RandomGenerator::runif(int size, int seed) {
    long double * arr = new long double[size];
    long double m = pow(2, 31) -1;
    arr[0] = seed;

    for(int i=0; i<=size;i++){
        arr[i+1] = RandomGenerator::LGMGenerator(m, arr[i]);
        arr[i] /= m;
    }
    return arr;
}

/* To generate Binomial Distribution
 *
 * @param:
 * int size: size of the binomial distribution generated
 * int n: n units wish to occur
 * double p: probability of the occurence in n units of uniform distribution: P(X <= p)
 * int seed: seed which is used to generate standard normal
 * */
long double* RandomGenerator::rbinom(int size, int n, double p, int seed) {
    int uniSize = size * n;

    long double * uniformArr = RandomGenerator::runif(uniSize, seed);
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


/* To generate Exponential Distribution
 *
 * @param:
 * long double *arr: array size
 * int size: size of array
 * double lambda: lambda in the expoenetial distribution
 * */
long double* RandomGenerator::rexp(long double *arr, int size, double lambda){

    long double *expArr = new long double[size];

    for(int i=0; i<size; i++){
        expArr[i] = -lambda * log(arr[i]);
    }

    return expArr;
}


/* Simulate Standard Normal Distribution Using Box-Muller
 *
 * @param:
 * long double *arr: uniform distribution array
 * int size: size of array
 * */
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


/* Simulate Standard Normal Distribution Using Polar-Marsaglia
 *
 * @param:
 * long double *arr: uniform distribution array
 * int size: size of array
 * int normSize: size of the normal distribution wish to be generated , normSize has to be
 * smaller than uniform distribution
 *
 * NOTE: For simplicity, the generated number of normal will be only half of the input uniform distribution
 * */
long double* RandomGenerator::polarmarsaglia(long double *arr, int size){

    int normSize = size /2;
    long double* normArr = new long double[normSize];
//    int arrSize = 1;
    int arrSize =0;
    long double v1 = 0;
    long double v2 = 0;
    long double w = 0;

//
//    for(int i=0; i<size; i+=2){
//        v1 = 2 * arr[i] - 1;
//        v2 = 2 * arr[i+1] -1;
//        w = v1 * v1  + v2 * v2; // using pow() would significantly increase the execution time.
//
//        if(w <= 1.0 ){
//            normArr[arrSize] = v1 * sqrt(-2 * log(w) / w);
//            normArr[arrSize+1] = v2 * sqrt(-2 * log(w) / w);
//            arrSize +=2;
//        }
//    }
//    //Put the total number of random variables generated in the first cell.
//    normArr[0] = arrSize -1;


    for(int i=0; i < normSize; i++){
        v1 = 2 * arr[i] - 1;
        v2 = 2 * arr[i+1] -1;
        w = v1 * v1  + v2 * v2; // using pow() would significantly increase the execution time.

        if(w <= 1.0 ){
            normArr[arrSize] = v1 * sqrt(-2 * log(w) / w);
            normArr[arrSize+1] = v2 * sqrt(-2 * log(w) / w);
            arrSize +=2;
        }
    }

    return normArr;
}


/* Generate the X of the bivariate normal
 *
 * @param:
 * long double *z1: standard normal array
 * double size: size of the array
 * */
long double* RandomGenerator::bivariateNormalX(long double *z1, int size){

    long double * arrX = new long double[size];

    double muX = 0;
    double sigmaX = 1;

    for(int i = 0; i<size; i++){
        arrX[i] = muX + sigmaX * z1[i];
    }

    return arrX;
}


/* Generate the Y in the bivariate normal
 *
 * @param:
 * long double *z1: standard normal array
 * long double *z2: standard normal array
 * double rho: correlation between z1 and z2
 * double size:  size of both z1 and z2
 * */
long double* RandomGenerator::bivariateNormalY(long double *z1, long double *z2, double rho, int size) {

    long double * arrY = new long double[size];

    double muY = 0;
    double sigmaY = 1;

    for(int i = 0; i<size; i++){
        arrY[i] = muY + sigmaY * rho * z1[i] + sigmaY * sqrt(1- rho * rho) * z2[i];
    }

    return arrY;
}

/* To generate Wiener process
 *
 * @param:
 * double t: Time period
 * int size: array size of t
 * int seed: seed which is used to generate standard normal
 * */
long double* RandomGenerator::wienerProcess(double t, int size, int seed){

    long double * stdNor= RandomGenerator::boxmuller(RandomGenerator::runif(size, seed), size);

    return Mutils::MatrixMultiply(stdNor, size, sqrt(t));
}
