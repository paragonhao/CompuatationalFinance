//
// Created by paragonhao on 22/1/19.
//

#include "OptionPricing.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
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

/* Call Option Price using Black-holes
 *
 * @param:
 * double r: interest rate
 * double sigma: sigma
 * double T: Period during which the process exists
 * double s0: initial stock price of x
 * double x: strike price
 * */
double OptionPricing::putOptionPriceBS(double r, double sigma, double T, double s0, double x){
    long double d1;
    long double d2;

    d1 = (log(s0/x) + (r+ 0.5*sigma*sigma) * T)/(sigma*sqrt(T));
    d2 = d1 - sigma * sqrt(T);

    return x *  Mutils::pnorm(-d2)/exp(r*T) - s0 * Mutils::pnorm(-d1);
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
double OptionPricing::callOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                       const double& sigma, const double& t, const int& steps){
    double delta = t/steps;
    double c = 0.0;
    double d = 0.0;
    double u = 0.0;
    double p_up = 0.0;
    double p_down = 0.0;
    double discount = exp(r*(delta));

    if(method == "a"){
        c = 0.5 * (exp(-r * delta) + exp((r + sigma * sigma)* delta));

        // volatility determines the up and down factors
        d = c - sqrt(c*c - 1);
        u = 1/d;
        p_up = (exp(r * delta) - d)/(u - d);
    }else if(method == "b"){
        u = exp(r * delta) * (1 + sqrt(exp(sigma * sigma * delta) - 1));
        d = exp(r * delta) * (1 - sqrt(exp(sigma * sigma * delta) - 1));
        p_up = 0.5;
    }else if(method == "c"){
        u = exp((r - sigma * sigma * 0.5) * delta + sigma * sqrt(delta));
        d = exp((r - sigma * sigma * 0.5) * delta - sigma * sqrt(delta));
        p_up = 0.5;
    }else if(method == "d"){
        u = exp(sigma * sqrt(delta));
        d = exp(-sigma * sqrt(delta));
        p_up = 0.5 + 0.5 * ((r - 0.5 * sigma * sigma) * sqrt(delta) / sigma);
    }

    p_down = 1 - p_up;

    vector< vector< double > > tree ( steps+1, vector<double> ( steps +1, 0 ) );
    vector< vector< double > > optionPriceTree ( steps+1, vector<double> ( steps +1, 0 ) );


    // generate the binomial tree based on the p_up, d, u factors
    for(int i=0; i<=steps; i++){
        for(int j=0; j<=i; j++){
            tree[i][j] = S * pow(u, j) * pow(d, i-j);
        }
    }

    for(int i=0; i<=steps; i++){
        for(int j=0; j<=i; j++){
            optionPriceTree[i][j] = max(0.0, tree[i][j] - K);
        }
    }

    for(int i = steps - 1; i >= 0; i--){
        for(int j = 0; j<= i; j++){
            optionPriceTree[i][j] = (p_down * optionPriceTree[i+1][j] + p_up * optionPriceTree[i+1][j+1])/discount;
        }
    }

    return optionPriceTree[0][0];

}

/* Call Option Pricing by using American binomial trees
 *
 * @param:
 * double S: spot price
 * double K: Strike price
 * double r: Interest rate
 * double Sigma: volatility
 * double t: time to maturity
 * double step: Number of periods for binomial trees
 * */
double OptionPricing::putOptionAmericanBinomial(const string& method, const double& S, const double& K, const double& r,
                                                 const double& sigma, const double& t, const int& steps){
    double delta = t/steps;
    double c = 0.0;
    double d = 0.0;
    double u = 0.0;
    double p_up = 0.0;
    double p_down = 0.0;
    double discount = exp(r*(delta));

    if(method == "a"){
        c = 0.5 * (exp(-r * delta) + exp((r + sigma * sigma)* delta));

        // volatility determines the up and down factors
        d = c - sqrt(c*c - 1);
        u = 1/d;
        p_up = (exp(r * delta) - d)/(u - d);
    }else if(method == "b"){
        u = exp(r * delta) * (1 + sqrt(exp(sigma * sigma * delta) - 1));
        d = exp(r * delta) * (1 - sqrt(exp(sigma * sigma * delta) - 1));
        p_up = 0.5;
    }else if(method == "c"){
        u = exp((r - sigma * sigma * 0.5) * delta + sigma * sqrt(delta));
        d = exp((r - sigma * sigma * 0.5) * delta - sigma * sqrt(delta));
        p_up = 0.5;
    }else if(method == "d"){
        u = exp(sigma * sqrt(delta));
        d = exp(-sigma * sqrt(delta));
        p_up = 0.5 + 0.5 * ((r - 0.5 * sigma * sigma) * sqrt(delta) / sigma);
    }

    p_down = 1 - p_up;

    vector< vector< double > > tree ( steps+1, vector<double> ( steps +1, 0 ) );
    vector< vector< double > > optionPriceTree ( steps+1, vector<double> ( steps +1, 0 ) );


    // generate the binomial tree based on the p_up, d, u factors
    for(int i=0; i<=steps; i++){
        for(int j=0; j<=i; j++){
            tree[i][j] = S * pow(u, j) * pow(d, i-j);
        }
    }

    for(int i=0; i<=steps; i++){
        for(int j=0; j<=i; j++){
            optionPriceTree[i][j] = max(0.0, K - tree[i][j]);
        }
    }

    for(int i = steps - 1; i >= 0; i--){
        for(int j = 0; j<= i; j++){
            optionPriceTree[i][j] = max((p_down * optionPriceTree[i+1][j] + p_up * optionPriceTree[i+1][j+1]),optionPriceTree[i][j])/discount;
        }
    }

    // print out the entire binomial tree
//    for(int i=0; i<=steps; i++){
//        for(int j=0; j<=i; j++){
//            cout<< tree[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//    for(int i=0; i<=steps; i++){
//        for(int j=0; j<=i; j++){
//            cout<< optionPriceTree[i][j] << " ";
//        }
//        cout << endl;
//    }

    return optionPriceTree[0][0];

}



/* Put Option Pricing by using binomial trees
 *
 * @param:
 * double S: spot price
 * double K: Strike price
 * double r: Interest rate
 * double Sigma: volatility
 * double t: time to maturity
 * double step: Number of periods for binomial trees
 * */
double OptionPricing::putOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                                 const double& sigma, const double& t, const int& steps){
    double delta = t/steps;
    double c = 0.0;
    double d = 0.0;
    double u = 0.0;
    double p_up = 0.0;
    double p_down = 0.0;
    double discount = exp(r*(delta));

    if(method == "a"){
        c = 0.5 * (exp(-r * delta) + exp((r + sigma * sigma)* delta));

        // volatility determines the up and down factors
        d = c - sqrt(c*c - 1);
        u = 1/d;
        p_up = (exp(r * delta) - d)/(u - d);
    }else if(method == "b"){
        u = exp(r * delta) * (1 + sqrt(exp(sigma * sigma * delta) - 1));
        d = exp(r * delta) * (1 - sqrt(exp(sigma * sigma * delta) - 1));
        p_up = 0.5;
    }else if(method == "c"){
        u = exp((r - sigma * sigma * 0.5) * delta + sigma * sqrt(delta));
        d = exp((r - sigma * sigma * 0.5) * delta - sigma * sqrt(delta));
        p_up = 0.5;
    }else if(method == "d"){
        u = exp(sigma * sqrt(delta));
        d = exp(-sigma * sqrt(delta));
        p_up = 0.5 + 0.5 * ((r - 0.5 * sigma * sigma) * sqrt(delta) / sigma);
    }

    p_down = 1 - p_up;

    vector< vector< double > > tree ( steps+1, vector<double> ( steps +1, 0 ) );
    vector< vector< double > > optionPriceTree ( steps+1, vector<double> ( steps +1, 0 ) );


    // generate the binomial tree based on the p_up, d, u factors
    for(int i=0; i<=steps; i++){
        for(int j=0; j<=i; j++){
            tree[i][j] = S * pow(u, j) * pow(d, i-j);
        }
    }

    for(int i=0; i<=steps; i++){
        for(int j=0; j<=i; j++){
            optionPriceTree[i][j] = max(0.0, K - tree[i][j]);
        }
    }

    for(int i = steps - 1; i >= 0; i--){
        for(int j = 0; j<= i; j++){
            optionPriceTree[i][j] = (p_down * optionPriceTree[i+1][j] + p_up * optionPriceTree[i+1][j+1])/discount;
        }
    }


    // print out the entire binomial tree
//    for(int i=0; i<=steps; i++){
//        for(int j=0; j<=i; j++){
//            cout<< tree[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//    for(int i=0; i<=steps; i++){
//        for(int j=0; j<=i; j++){
//            cout<< optionPriceTree[i][j] << " ";
//        }
//        cout << endl;
//    }


    return optionPriceTree[0][0];

}