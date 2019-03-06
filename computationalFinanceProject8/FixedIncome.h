//
// Created by paragonhao on 3/3/19.
//

#ifndef COMPUTATIONALFINANCEPROJECT8_FIXEDINCOME_H
#define COMPUTATIONALFINANCEPROJECT8_FIXEDINCOME_H
#include <Eigen/Dense>
#include <ctime>

using namespace Eigen;
using namespace std;

class FixedIncome {

public:
    static double calculateZCBPV(double r0, double sigma, double specialK, double rbar, double FV, double T, int steps);

    static double getBondPrice(double FV, double rt, double specialK, double sigma, double rbar, double T, double t);

    static void getRPathAndBigRVasicek(double r0, double sigma, double specialK, double rbar, double T, int steps,
                                       MatrixXd &rMat, VectorXd &bigR, int simNum);

    static void getRPathAndBigRCIR(double r0, double sigma, double specialK, double rbar, double T, int steps, MatrixXd &rMat,
                            VectorXd &bigR, int simNum);

    static void getfunctionAandB(double specialK, double sigma, double rbar, double t, double T, double &A_t_T, double &B_t_T);

    static void getRPathAndBigRGPlusPlusModel(int simNum, double T, int steps, double a, double b, double sigma, double phi, double eta,
                                              double r0, double rho, MatrixXd &rMat, VectorXd & bigR, MatrixXd &xMat,MatrixXd &yMat);

    static void
    getExplicitOptionCallPriceCIR(double r0, double sigma, double specialK, double rbar, double strike, double FV,
                                  double T, double S, double t);
};


#endif //COMPUTATIONALFINANCEPROJECT8_FIXEDINCOME_H
