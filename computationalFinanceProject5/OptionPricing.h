//
// Created by paragonhao on 22/1/19.
//

#ifndef COMPUTATIONALFINANCE_OPTIONPRICING_H
#define COMPUTATIONALFINANCE_OPTIONPRICING_H
#include <Eigen/Dense>
#include <string>

using namespace std;
using namespace Eigen;

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
    static double stockPriceAtT(const double& S, const double& r, const double& sigma, const double& t, const double &wt);

    static double putOptionEuropeanBinomial(const string& method, const double& S, const double& K, const double& r,
                                             const double& sigma, const double& t, const int& steps);
    static double putOptionPriceBS(double r, double sigma, double T, double s0, double x);
    static double putOptionAmericanBinomial(const string& method, const double& S, const double& K, const double& r,
                                            const double& sigma, const double& t, const int& steps);
    static void generatePricePath(const double *tArray,  const double & nSims, const double & halfPath, const double &s, const double &r,
                                  const double &sigma, ArrayXXd &priceProcess);
    static void calculateMatrixA(ArrayXXd &priceProcess, MatrixXd &matA,
                                 const int &k, const int &nPath, const int &currentSimCol, const int &method);
    static void calcualateMatrixb(ArrayXXd & cashflowMatrix, ArrayXXd & priceProcess, VectorXd & b, ArrayXXd & index,
                                  const int &k, const int &nPath, const int &currentSimCol, const int &method, const int &nSims,
                                  const double &r, const double &delta, const double x);
    static void calculateContinuationValue(VectorXd &a, ArrayXXd & priceProcess, ArrayXXd & index,
                                           ArrayXXd & cashFlowMatrix, const int &nSims,
                                           const int &currentSimCol, const int &nPath, const int &k, const int &method);
    static double calculateFinalPayOff(ArrayXXd & cashFlowMatrix, ArrayXXd &index,
                                       const int &nPath, const int &nSims, const double &r, const double & delta);
};


#endif //COMPUTATIONALFINANCE_OPTIONPRICING_H
