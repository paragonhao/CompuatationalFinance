#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <fstream>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"

using namespace std;

inline double alphaX(double x_t){
    return 0.2 - 0.5 * x_t;
}

inline double betaX(double x_t, double w_t){
    return 2.0/3 * w_t;
}

inline double alphaY(double delta_t, double y){
    return (2.0/(1.0 + delta_t)) * y  + ((1 + delta_t * delta_t * delta_t) / 3.0);
}

inline double solveQn1DXT(double t, int size, int seed){
    double delta_t = t /size;
    double x = 1;

    long double * dW_t = RandomGenerator::wienerProcess(delta_t, size, seed);

    for (int i =0 ; i< size ; i++){
        x = x + alphaX(x) * delta_t + betaX(x, dW_t[i]);
    }
    return x;
}

inline double solveQn1DYT(double t, int size, int seed){
    double delta_t = t/size;
    double y = 0.75;

    long double * dZ_t = RandomGenerator::wienerProcess(delta_t, size, seed);

    for (int i = 0 ; i < size; i++){
        double curr_delta = delta_t * i;
        y = y + double(alphaY(curr_delta, y) * delta_t + ((1 + curr_delta * curr_delta * curr_delta) / 3.0) * dZ_t[i]);
    }
    return y;
}

inline double solveQn2DXT(double t, int size, int seed){
    double delta_t = t/size;
    double x = 1.0;

    long double * dW_t = RandomGenerator::wienerProcess(delta_t, size*2, seed);
    long double * dZ_t = &dW_t[size];

    for(int i = 0; i<size; i++){
        // Euler
        // x = x +  double(0.25 * x * delta_t + (1.0/3) * x * dW_t[i] - 0.75 * x * dZ_t[i]);
        // Milshtein's Scheme
        x = x +  double(0.25 * x * delta_t + (1.0/3) * x * dW_t[i] - 0.75 * x * dZ_t[i]
                + 0.5 * (1.0/3) * x * (1.0/3) * (dW_t[i] * dW_t[i] - delta_t)
                        + 0.5 * 0.75 * 0.75 * x * (dZ_t[i] * dZ_t[i] - delta_t));
    }
    return x;

}

inline long double * solveQn2YT(double t, int size, int seed){

    long double * w_t = RandomGenerator::wienerProcess(t, size * 2, seed);
    long double * z_t = &w_t[size];
    long double * y_t = new long double [size];

    for(int i=0; i < size; i++){
        y_t[i] = cbrt(exp(-0.08 * t + (1.0/3) * w_t[i] + 0.75 * z_t[i]) + 1);
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
    cout << "E[X_2^(1/3)]: " << x1_3/size << endl;


    // Find Y3_T based on dY_t = ((2/1 + t) * Y_t + (1+t^3) / 3) * dt + (1 + t^3)/3 * dZ_t
    long double y3Arr[size];
    t = 3;

    for(int i = 0; i < size ; i ++){
        y3Arr[i] = solveQn1DYT(t, size, seed + i * 100);
    }
    cout << "E[Y_3]: " << Mutils::Mean(y3Arr, size) << endl;


    // Find Y2_T based on dY_t = ((2/1 + t) * Y_t + (1+t^3) / 3) * dt + (1 + t^3)/3 * dZ_t
    long double y2Arr[size];
    t = 2;

    for(int i = 0; i < size ; i++){
        y2Arr[i] = solveQn1DYT(t, size, seed + i * 33);
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
    int size = 1000;
    long double x_t[size];
    long double *y_t;
    double t = 3;

    for(int i= 0; i < size; i++) {
        x_t[i] = cbrt(solveQn2DXT(t, size, seed + i * 100) + 1.0);
    }

    cout << "E[(1+X_3)^1/3] is : "<< Mutils::Mean(x_t, size) << endl;


    size = 1000;
    y_t = solveQn2YT(t, size, seed);
    cout << "E[(1+Y_3)^1/3] is : "<< Mutils::Mean(y_t, size)<< endl;

}

void RunQn3(int seed){
    double r = 0.04;
    double sigma = 0.2;
    double s0 = 88;
    double T = 5;
    int size = 1000;
    double x = 100;


    cout << "3 (a) Function <callOptionSimAntithetic> is available at OptionPricing.cpp "<<endl;
    cout << "Assuming r = 0.04, sigma =0.2, s0= 88, T = 5, x= 100 for 1000 simulation"<<endl;
    long double * price = OptionPricing::callOptionSimAntithetic(seed,size,r,sigma,T,s0, x);
    cout <<  "Price using Antithetic Variates in Monte Carlo Simulation: $"<< Mutils::Mean(price, size) << endl;
    cout << endl;
    cout << "3 (b) Function <callOptionPriceBS> is available at OptionPricing.cpp "<<endl;
    double bsprice = OptionPricing::callOptionPriceBS(r, sigma, T, s0, x);
    cout <<  "Price using black scholes formula: $"<< bsprice << endl;
    cout <<endl;
    cout << "#################################### Qn 3(c) ###################################"<<endl;
    // TODO: TO be determined whether to use formula or approximation
}

void RunQn4(int seed){
    int size = 1000;
    auto rho = 0.6;
    auto r =0.03;
    auto s_0 = 48;
    auto v_0 = 0.05;
    auto sigma = 0.42;
    auto alpha = 5.8;
    auto beta = 0.0625;

    long double * z1;
    long double * z2;

    long double * x;
    long double * y;

    //generate normal distribution
    z1 = RandomGenerator::boxmuller(RandomGenerator::runif(size * 2, seed), size * 2);
    z2 = &z1[size];

    x = RandomGenerator::bivariateNormalX(z1, size);
    y = RandomGenerator::bivariateNormalY(z1, z2, rho, size);



}

int main() {
    int seed = 1234567890;
    cout << "#################################### Qn 1 ###################################" << endl;
    RunQn1(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 2 ###################################" << endl;
    RunQn2(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 3 ###################################" << endl;
    RunQn3(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 4 ###################################" << endl;
    RunQn4(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    return 0;
}

