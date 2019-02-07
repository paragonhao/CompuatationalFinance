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

/* Put Option Pricing by using trinomial trees
 *
 * @param:
 * string method: method chose to generate p_up, p_down, u & d
 * double S: spot price
 * double K: Strike price
 * double r: Interest rate
 * double Sigma: volatility
 * double t: time to maturity
 * double step: Number of periods for binomial trees
 * */
double OptionPricing::callOptionEuropeanTrinomial(const string& method, const double& S, const double& K, const double& r,
                                   const double& sigma, const double& t, const int& steps){
    double delta = t/steps;
    double d = 0.0;
    double u = 0.0;
    double p_up = 0.0;
    double p_down = 0.0;
    double p_m = 0.0;
    double discount = exp(r*(delta));
    double delta_xu = 0.0;
    double delta_xd = 0.0;

    if(method == "a"){
        d = exp(-sigma * sqrt(3 * delta));
        u = 1/d;
        p_down = (r*delta*(1-u) + (r*delta) *(r*delta) + sigma*sigma*delta)/((u-d)*(1-d));
        p_up = (r*delta*(1-d) + r*delta*r*delta + sigma*sigma*delta)/((u - d) *(u- 1));
    }else if(method == "b"){
        delta_xu = sigma * sqrt(3*delta);
        delta_xd = -delta_xu;
        double tempNumerator1 = sigma * sigma * delta + pow((r - 0.5 * sigma * sigma), 2)  * delta * delta;
        double tempNumerator2 = (r - 0.5 * sigma * sigma) * delta;
        p_down = 0.5*((tempNumerator1/(delta_xu*delta_xu)) - (tempNumerator2/delta_xu));
        p_up = 0.5*((tempNumerator1/(delta_xu*delta_xu)) + (tempNumerator2/delta_xu));
    }
    p_m = 1 - p_down - p_up;

    vector< vector< double > > tree ( steps+1, vector<double> ( 2 * steps +1, 0 ) );
    vector< vector< double > > optionPriceTree ( steps+1, vector<double> ( 2 * steps +1, 0 ) );

    if(method == "a"){
        // generate the trinomial tree based on the u, d factors
        for(int i=0; i<=steps; i++){
            for(int j=0; j<2*i + 1; j++){
                tree[i][j] = S * pow(u, max(i-j, 0)) * pow(d, max(j-i, 0));
            }
        }
        for(int i=0; i<=steps; i++){
            for(int j=0; j<2*i + 1; j++){
                optionPriceTree[i][j] = max(0.0, tree[i][j] - K);
            }
        }

    }else if(method == "b"){
        double X = log(S);
        for(int i=0; i<=steps; i++){
            for(int j=0; j<2*i + 1; j++){
                tree[i][j] = X  + delta_xu * max(i-j, 0) + delta_xd * max(j-i, 0);
            }
        }
        for(int i=0; i<=steps; i++){
            for(int j=0; j<2*i + 1; j++){
                optionPriceTree[i][j] = max(0.0, exp(tree[i][j]) - K);
            }
        }
    }




    for(int i = steps - 1; i >= 0; i--){
        int num = (i+1)*2 + 1;
        for(int j = 0; j< num; j++){
            if(j+2 <num){
                optionPriceTree[i][j] = (p_up * optionPriceTree[i+1][j]
                                         + p_m * optionPriceTree[i+1][j+1] + p_down * optionPriceTree[i+1][j+2])/discount;
            }else{
                break;
            }
        }
    }


//  print out the tree
//    for(int i=0; i<=steps; i++){
//        for(int j=0; j<2*i + 1; j++){
//            cout<< optionPriceTree[i][j]<<", ";
//        }
//        cout << endl;
//    }
//
//


    return optionPriceTree[0][0];
}


double OptionPricing::callOptionEuropeanLDS(const double& S, const double& K, const double& r,
                                                  const double& sigma, const double& T, const int& N,
                                            const int& base1, const int& base2){
    long double * base1Arr = RandomGenerator::getHaltonSequence(base1, N);
    long double * base2Arr = RandomGenerator::getHaltonSequence(base2, N);
    long double * NormArr;

    // Generate standardard normal using halton sequence
    NormArr = RandomGenerator::boxmullerWithHaltonSeq(base1Arr, base2Arr, N);
    // Generate wiener process
    long double * pay_off = OptionPricing::callOptionPriceSimulation(Mutils::MatrixMultiply(NormArr,N, sqrt(T)),N,r,sigma,T,S,K);

    return Mutils::Mean(pay_off, N);
}