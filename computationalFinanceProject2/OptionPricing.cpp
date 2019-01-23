//
// Created by paragonhao on 22/1/19.
//

#include "OptionPricing.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include "Mutils.h"


using namespace std;

long double *  OptionPricing::callOptionPriceSimulation(long double *w_t, int size, double r, double sigma, double T, double s0, double x){

    long double * s_T = new long double[size];

    for(int i = 0; i< size; i++){
        double price = double(s0 * exp((r - (sigma * sigma)/2) * T + sigma * w_t[i]) - x);

        s_T[i] = max(price, 0.0)/(exp(r*T));
    }
    return s_T;
}

double OptionPricing::callOptionPriceBS(double r, double sigma, double T, double s0, double x){
    double d1;
    double d2;

    d1 = (log(s0/x) + (r+ 0.5*sigma*sigma) * T)/(sigma*sqrt(T));
    d2 = d1 - sigma * sqrt(T);

    return s0* Mutils::pnorm(d1) - (x *  Mutils::pnorm(d2)/exp(r*T));
}
