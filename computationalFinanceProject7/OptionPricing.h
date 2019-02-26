//
// Created by paragonhao on 22/1/19.
//

#ifndef COMPUTATIONALFINANCE_OPTIONPRICING_H
#define COMPUTATIONALFINANCE_OPTIONPRICING_H

#include <string>

using namespace std;

class OptionPricing {
public:
    static double * callOptionPriceSimulation(double *w_t, int size, double r, double sigma, double T, double s0, double x);
    static double callOptionSimAntithetic(int seed, int size, double r, double sigma, double T, double s0, double x);
    static double callOptionPriceBS(double r, double sigma, double T, double s0, double x);
    static double callOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                                  const double& sigma, const double& t, const int& steps);
    static double callOptionEuropeanTrinomial(const string& method, const double& S, const double& K, const double& r,
                                              const double& sigma, const double& t, const int& steps);

    static double callOptionEuropeanLDS(const double& S, const double& K, const double& r,
                                        const double& sigma, const double& T, const int& N,
                                        const int& base1, const int& base2);


    static double putOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                             const double& sigma, const double& t, const int& steps);
    static double putOptionPriceBS(double r, double sigma, double T, double s0, double x);
    static double putOptionAmericanBinomial(const string& method, const double& S, const double& K, const double& r,
                                            const double& sigma, const double& t, const int& steps);
    static double FixedStrikeLookBackCall(const double &r, const double &T, const double &s0, const double &strike, const double &sigma,
                                                  const int & num, const int & steps, const double & delta);
    static double FixedStrikeLookBackPut(const double &r, const double &T, const double &s0, const double &strike, const double &sigma,
                                                 const int & num, const int & steps, const double & delta);
};


#endif //COMPUTATIONALFINANCE_OPTIONPRICING_H
