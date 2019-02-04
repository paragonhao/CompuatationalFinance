//
// Created by paragonhao on 22/1/19.
//

#include "OptionPricing.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include "Mutils.h"
#include "RandomGenerator.h"

using namespace std;

/* Call Option Price Simulation
 *
 * @param:
 * long double* : Wiener process generated
 * int : size of the array
 * double r: interest rate
 * double sigma: sigma
 * double T: Period during which the process exists
 * double s0: initial stock price of x
 * double x: strike price
 * */
long double * OptionPricing::callOptionPriceSimulation(long double *w_t, int size, double r, double sigma, double T, double s0, double x){

    long double * s_T = new long double[size];

    for(int i = 0; i< size; i++){
        double price = double(s0 * exp((r - (sigma * sigma)/2) * T + sigma * w_t[i]) - x);

        s_T[i] = max(price, 0.0)/(exp(r*T));
    }
    return s_T;
}

/* Call Option Price using Black-holes
 *
 * @param:
 * double r: interest rate
 * double sigma: sigma
 * double T: Period during which the process exists
 * double s0: initial stock price of x
 * double x: strike price
 * */
double OptionPricing::callOptionPriceBS(double r, double sigma, double T, double s0, double x){
    long double d1;
    long double d2;

    d1 = (log(s0/x) + (r+ 0.5*sigma*sigma) * T)/(sigma*sqrt(T));
    d2 = d1 - sigma * sqrt(T);

    return s0 * Mutils::pnorm(d1) - (x *  Mutils::pnorm(d2)/exp(r*T));
}


/* Call Option Price Simulation by using Antithetic Variance Reduction Method
 *
 * @param:
 * int seed : seed
 * int size: size of the array
 * double r: interest rate
 * double sigma: sigma
 * double T: Period during which the process exists
 * double s0: initial stock price of x
 * double x: strike price
 * */
long double OptionPricing::callOptionSimAntithetic(int seed, int size, double r, double sigma, double T, double s0, double x){

    long double *w_t = RandomGenerator::wienerProcess(T,size,seed);

    long double *s_T = OptionPricing::callOptionPriceSimulation(w_t, size, r, sigma, T, s0, x);

    long double * y = new long double[size];

    std::copy(w_t, w_t + size, y);
    Mutils::MatrixMultiply(y, size, -1); // we know that the variance is 5 in this case and expectaion value is 0

    long double *s_T_2 = OptionPricing::callOptionPriceSimulation(y, size, r, sigma, T, s0, x);


    long double *tArr = new long double[size];

    for(int i=0; i<size; i++){
        tArr[i] = 0.5 * (s_T[i] + s_T_2[i]);
    }
    return Mutils::Mean(tArr, size);
}


/* Call Option Pricing by using binomial trees
 *
 * @param:
 * double S: spot price
 * double K: Strike price
 * double r: Interest rate
 * double Sigma: volatility
 * double t: time to maturity
 * double step: Number of periods for binomial trees
 * */
double callOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                       const double& sigma, const double& t, const int& steps){
    double delta = t/steps;
    double c = 0.0;


    if(method == "a"){
       cout <<endl;
    }



    return 0.0;

}