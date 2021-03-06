//
// Created by paragonhao on 6/3/19.
//

#ifndef COMPUTATIONALFINANCEPROJECT9_MORTGAGEBACKEDSECURITIES_H
#define COMPUTATIONALFINANCEPROJECT9_MORTGAGEBACKEDSECURITIES_H


#include <Eigen/Dense>

using namespace Eigen;

class MortgageBackedSecurities {

public:
    static double getBurnOutRate(const double &pv0, const double &pvt_minus_1);

    static double getInterestRate(const double &annualizedMortgageRate, const double &tenYearRateT_minus_1);

    static double getSeasoningSG(const double &t);

    static double getSeasonalitySY(const int &month_index);

    static void getRPathCIRModel(double r0, double sigma, double specialK, double rbar, double T, int steps,
                         MatrixXd &rMat, int simNum);

    static double getCPR_t(const double &pv0, const double &pvt_minus_1, const double &t, const double &R, const double &rt,
                    const int &monthIndex);

    static double getNumerixPrepaymentModel(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double duration, double OAS, bool isOAS);
};


#endif //COMPUTATIONALFINANCEPROJECT9_MORTGAGEBACKEDSECURITIES_H
