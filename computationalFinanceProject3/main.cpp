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

inline long double reflectionTruncationModel(int seed, double T, int size, int x, double rho, double r,
                                             double s0, double v0, double sigma, double a, double b){
    double delta = T/size;
    long double * z1;
    long double * z2;

    long double * w1;
    long double * w2;


    double v_reflect = v0;
    double s_reflect = s0;

    //generate normal distribution
    z1 = RandomGenerator::boxmuller(RandomGenerator::runif(size * 2, seed), size * 2);
    z2 = &z1[size];

    w1 = RandomGenerator::bivariateNormalX(z1, size);
    w2 = RandomGenerator::bivariateNormalY(z1, z2, rho, size);

    for(int i = 0; i < size; i++){
        s_reflect = s_reflect + r * s_reflect * delta + sqrt(abs(v_reflect)) * s_reflect * w1[i] * sqrt(delta);
        v_reflect = abs(v_reflect) + a * (b - abs(v_reflect)) * delta + sigma * sqrt(abs(v_reflect))* w2[i] * sqrt(delta);
    }
    return  max(s_reflect - x, 0.0);
}

inline double partialTruncationModel(int seed, double T, int size, int x, double rho, double r,
                                  double s0, double v0, double sigma, double a, double b){
    double delta = T/size;
    long double * z1;
    long double * z2;

    long double * w1;
    long double * w2;

    double v_partial = v0;
    double s_partial = s0;

    //generate normal distribution
    z1 = RandomGenerator::boxmuller(RandomGenerator::runif(size * 2, seed), size * 2);
    z2 = &z1[size];

    w1 = RandomGenerator::bivariateNormalX(z1, size);
    w2 = RandomGenerator::bivariateNormalY(z1, z2, rho, size);

    for(int i = 0; i < size; i++){
        s_partial += r * s_partial * delta + sqrt(max(0.0, v_partial)) * s_partial * w1[i] * sqrt(delta);
        v_partial += a * (b - v_partial) * delta + sigma * sqrt(max(0.0, v_partial)) * w2[i] * sqrt(delta);
    }

    return  max(s_partial - x, 0.0);
}

inline double fullTruncationModel(int seed, double T, int size, int x, double rho, double r,
                                double s0, double v0, double sigma, double a, double b){

    double delta = T/size;
    long double * z1;
    long double * z2;

    long double * w1;
    long double * w2;

    double v_full = v0;
    double s_full = s0;

    //generate normal distribution
    z1 = RandomGenerator::boxmuller(RandomGenerator::runif(size * 2, seed), size * 2);
    z2 = &z1[size];

    w1 = RandomGenerator::bivariateNormalX(z1, size);
    w2 = RandomGenerator::bivariateNormalY(z1, z2, rho, size);

    for(int i = 0; i < size; i++){
        s_full +=r * s_full * delta + sqrt(max(0.0, v_full)) * s_full * w1[i] * sqrt(delta);
        v_full += a * (b - max(0.0,v_full)) * delta + sigma * sqrt(max(0.0, v_full)) * w2[i] * sqrt(delta);
    }
    return  max(s_full - x, 0.0);
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

    int iter = 11;
    double x = 20;
    double sigma = 0.25;
    double r = 0.04;
    double T = 0.5;
    double epsilon = 0.01;
    long double greek_delta[iter];
    long double greek_theta[iter];
    long double greek_vega[iter];
    long double greek_rho[iter];
    long double greek_gamma[iter];

    int s0 = 15;
    int size = 1000;

    for(int i = 0; i < iter; i++) {

        double price = OptionPricing::callOptionSimAntithetic(seed,size,r,sigma,T,s0, x);

        greek_delta[i] = (OptionPricing::callOptionSimAntithetic(seed,size,r,sigma,T,s0 + epsilon, x) -
                price)/ epsilon;

        greek_theta[i] = (OptionPricing::callOptionSimAntithetic(seed,size,r,sigma,T + epsilon,s0, x) -
                price)/epsilon;

        greek_vega[i] = (OptionPricing::callOptionSimAntithetic(seed,size,r,sigma + epsilon * epsilon, T, s0, x) -
                price) / (epsilon * epsilon);

        greek_rho[i] = (OptionPricing::callOptionSimAntithetic(seed,size,r + epsilon * epsilon,sigma , T, s0, x) -
                price) / (epsilon * epsilon);

        greek_gamma[i] = (OptionPricing::callOptionSimAntithetic(seed, size, r, sigma , T, s0 + epsilon * 2, x) -
                2 * OptionPricing::callOptionSimAntithetic(seed, size, r, sigma , T, s0 + epsilon, x) + price)
                         / (epsilon*epsilon);

        s0++;
    }
    cout << "Delta: ";
    for(int i=0; i< iter;i++){
        cout <<greek_delta[i] << ", ";
    }
    cout << endl;

    cout << "Gamma: ";
    for(int i=0; i< iter;i++){
        cout <<greek_gamma[i] << ", ";
    }
    cout << endl;

    cout << "Theta: ";
    for(int i=0; i< iter;i++){
        cout <<greek_theta[i] << ", ";
    }
    cout << endl;

    cout << "Rho: ";
    for(int i=0; i< iter;i++){
        cout <<greek_rho[i] << ", ";
    }
    cout << endl;

    cout << "Vega: ";
    for(int i=0; i< iter;i++){
        cout <<greek_vega[i] << ", ";
    }
    cout << endl;


}

void RunQn4(int seed){
    int size = 1000;
    double rho = -0.6;
    double r = 0.03;
    double s_0 = 48;
    double v_0 = 0.05;
    double sigma = 0.42;
    double alpha = 5.8;
    double beta = 0.0625;
    double T = 0.5;
    double x = 50;

    long double * full = new long double[size];
    long double * partial = new long double[size];
    long double * reflection = new long double[size];

    for(int i =0; i< size; i++){

        full[i] = fullTruncationModel(seed + i * 1000, T, size, x, rho, r, s_0, v_0, sigma, alpha, beta);
        partial[i] = partialTruncationModel(seed + i * 1000, T, size, x, rho, r, s_0, v_0, sigma, alpha, beta);
        reflection[i] = reflectionTruncationModel(seed + i * 1000, T, size, x, rho, r, s_0, v_0, sigma, alpha, beta);
    }
    cout << Mutils::Mean(full,size) << endl;
    cout << Mutils::Mean(partial,size) << endl;
    cout << Mutils::Mean(reflection,size) << endl;
}

int main() {
    int seed = 1234567890;
    cout << "#################################### Qn 1 ###################################" << endl;
    //RunQn1(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 2 ###################################" << endl;
    //RunQn2(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 3 ###################################" << endl;
    //RunQn3(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 4 ###################################" << endl;
    RunQn4(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    return 0;
}

