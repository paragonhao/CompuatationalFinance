#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <fstream>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"

using namespace std;

double alphaX(double x_t){
    return 0.2 - 0.5 * x_t;
}

double betaX(double x_t, double w_t){
    return 2.0/3 * w_t;
}

double solveQn1DXT(double t, int size, int seed){
    double delta_t = t /size;
    double x = 1;

    long double * dW_t = RandomGenerator::wienerProcess(delta_t, size, seed);

    for (int i =0 ; i< size ; i++){
        x = x + alphaX(x) * delta_t + betaX(x, dW_t[i]);
    }
    return x;
}

double alphaY(double delta_t, double y){
    return 2.0/(1.0 + delta_t) * y  + (1 + delta_t * delta_t * delta_t) / 3.0;
}

double solveQn1DYT(double t, int size, int seed){
    double delta_t = t/size;
    double y = 0.75;

    long double * dZ_t = RandomGenerator::wienerProcess(delta_t, size, seed);

    for (int i = 0 ; i < size; i++){
        y = y + alphaY(delta_t * (i + 1), y) * delta_t + ((1 + delta_t * delta_t * delta_t) / 3.0) * dZ_t[i];
    }
    return y;
}

double solveQn2DXT(double t, int size, int seed){
    double delta_t = t/size;
    double x = 1.0;

    long double * dW_t = RandomGenerator::wienerProcess(delta_t, size, seed);
    long double * dZ_t = RandomGenerator::wienerProcess(delta_t, size, seed + 1000);

    for(int i = 0; i<size ; i++){
        x = x +  double(0.25 * x * delta_t +  (1.0/3) * x * dW_t[i] - 0.75 * x * dZ_t[i]
        + 0.5 * (1.0/3) * x * (1.0/3) * (dW_t[i] * dW_t[i] - delta_t)
                + 0.5 * 0.75 * 0.75 * x * (dZ_t[i] * dZ_t[i] - delta_t));
    }
    return x;

}

long double * solveQn2YT(double t, int size, int seed){

    long double * w_t = RandomGenerator::wienerProcess(t, size, seed);
    long double * z_t = RandomGenerator::wienerProcess(t, size, seed + 1000);
    long double * y_t = new long double [size];

    for(int i=0; i < size; i++){
        y_t[i] = cbrt(exp(-0.008 * t + (1.0/3) * w_t[i] + 0.75 * z_t[i]) + 1);
    }

    return y_t;
}

// Using Euler scheme
void RunQn1(int seed){
    int size = 2000;

    //  Find X_T based on dX_t = (1/5 - 1/2 X_t) * dt + 2/3 * dW_t
    long double * x2Arr = new long double[size];
    double t = 2;

    double x1_3 = 0;
    for(int i= 0; i < size; i++) {
        x2Arr[i] = solveQn1DXT(t, size, seed + i);
        x1_3 += cbrt(x2Arr[i]);
    }
    cout << "E[X_2^(1/3)]: " << x1_3/2000  << endl;


    // Find Y3_T based on dY_t = ((2/1 + t) * Y_t + (1+t^3) / 3) * dt + (1 + t^3)/3 * dZ_t
    long double y3Arr[size];
    t = 3;

    for(int i = 0; i < size ; i ++){
        y3Arr[i] = solveQn1DYT(t, size, seed + i);
    }
    cout << "E[Y_3]: " << Mutils::Mean(y3Arr, size) << endl;


    // Find Y2_T based on dY_t = ((2/1 + t) * Y_t + (1+t^3) / 3) * dt + (1 + t^3)/3 * dZ_t
    long double y2Arr[size];
    t = 2;

    for(int i = 0; i < size ; i++){
        y2Arr[i] = solveQn1DYT(t, size, seed + i);
    }

    int counter = 0;
    for(int i =0; i < size; i++){
        double temp = (y2Arr[i] > 5) ? 1: 0;
        counter += temp;
    }
    cout << "P(Y_2 > 5): " << double(counter)/size<< endl;

    // Find E[X_2* Y_2 * 1(X_2 > 1)]
    // we already have computed X_2 and Y_2
    long double ex2y2 = 0;
    for(int i=0; i < size; i++){
        double temp = (x2Arr[i] >1)? 1 : 0;
        ex2y2 += x2Arr[i] * y2Arr[i] * temp;
    }
    cout << "E[X_2Y_21(X_2 > 1)] : " << ex2y2/size << endl;
}


// Using Milshtein scheme
void RunQn2(int seed){
    int size = 2000;
    long double x_t[size];
    long double *y_t;
    double t = 3;


    for(int i= 0; i < size; i++) {
        x_t[i] = cbrt(solveQn2DXT(t, size, seed + i) + 1.0);
    }

    double x_mean = Mutils::Mean(x_t, size);
    cout << "E[(1+X_3)^1/3] is : "<< x_mean << endl;

    y_t = solveQn2YT(t, size, seed);
    cout << "E[(1+Y_3)^1/3] is : "<<Mutils::Mean(y_t, size)<< endl;

}



int main() {
    int seed = 1234567890;
    //RunQn1(seed);
    RunQn2(seed);
    return 0;
}

