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
    static double getInterestRateSim(double r0, double sigma, double specialK, double rbar, double T, int steps, double payoff);

    static double getBondPrice(double FV, double rt, double specialK, double sigma, double rbar, double T, double t);
    static void getRPathAndBigR(double r0, double sigma, double specialK, double rbar, double T, int steps, MatrixXd & rMat, VectorXd & bigR);
};


#endif //COMPUTATIONALFINANCEPROJECT8_FIXEDINCOME_H
