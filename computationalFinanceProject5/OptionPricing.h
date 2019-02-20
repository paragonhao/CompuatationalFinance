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
    static void generatePricePath(const double * tArray,  const double & nSims, const double & halfPath, const double &s, const double &r,
                                    const double &sigma, vector< vector< double > >  &priceProcess);
    static void calculateMatrixA(vector< vector< double > > &priceProcess, MatrixXd &matA,
                                 const int &k, const int &nPath, const int &currentSimCol, const int &method);
    static void calcualateMatrixb(vector< vector< double > > & cashflowMatrix, vector< vector< double > > & priceProcess, VectorXd & b,vector< vector< double > > & index,
                                  const int &k, const int &nPath, const int &currentSimCol, const int &method, const int &nSims,
                                  const double &r, const double &delta, const double x);
    static void calculateContinuationValue(VectorXd &a, vector< vector< double > > & priceProcess, vector< vector< double > > & index,
                                                   vector< vector< double > > & cashFlowMatrix, const int &nSims,
                                                   const int &currentSimCol, const int &nPath, const int &k, const int &method);
    static double calculateFinalPayOff(vector< vector< double > > & cashFlowMatrix, vector< vector< double > > &index,
                                       const int &nPath, const int &nSims, const double &r, const double & delta);
};


#endif //COMPUTATIONALFINANCE_OPTIONPRICING_H
