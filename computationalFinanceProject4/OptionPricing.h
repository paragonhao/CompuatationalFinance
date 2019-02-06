//
// Created by paragonhao on 22/1/19.
//

#ifndef COMPUTATIONALFINANCE_OPTIONPRICING_H
#define COMPUTATIONALFINANCE_OPTIONPRICING_H

#include <string>

using namespace std;

class OptionPricing {
public:
    static long double * callOptionPriceSimulation(long double *w_t, int size, double r, double sigma, double T, double s0, double x);
    static long double callOptionSimAntithetic(int seed, int size, double r, double sigma, double T, double s0, double x);
    static double callOptionPriceBS(double r, double sigma, double T, double s0, double x);
    static double callOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                                  const double& sigma, const double& t, const int& steps);
    static double callOptionEuropeanTrinomial(const string& method, const double& S, const double& K, const double& r,
                                              const double& sigma, const double& t, const int& steps);

    static double putOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                             const double& sigma, const double& t, const int& steps);
    static double putOptionPriceBS(double r, double sigma, double T, double s0, double x);
    static double putOptionAmericanBinomial(const string& method, const double& S, const double& K, const double& r,
                                            const double& sigma, const double& t, const int& steps);
};


#endif //COMPUTATIONALFINANCE_OPTIONPRICING_H
