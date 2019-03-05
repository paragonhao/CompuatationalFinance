#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"
#include "JDDefaultOption.h"
#include "DifferenceMethod.h"
#include "FixedIncome.h"
#include <random>
#include <Eigen/Dense>
#include "include/stats.hpp"

using namespace Eigen;
using namespace std;

void RunQn1a(double r0, double sigma, double specialK, double rbar){
    double FV = 1000;
    double T = 0.5;
    int steps = 127;
    cout << "The bond price at time 0 is: "<<  FixedIncome::calculateZCBPV(r0, sigma, specialK, rbar, FV, T, steps) << endl;
}

void RunQn1b(double r0, double sigma, double specialK, double rbar){
    int numCF =8;
    vector<double> cashflow{ 30, 30, 30, 30, 30, 30, 30, 1030};
    vector<double> time{0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4};
    int tradingDays = 252;

    // treat each dividend payment as a zero coupon bond,
    double pv = 0;
    for(int i=0; i< numCF; i++){
        auto steps = floor(time[i] * tradingDays);
        pv += FixedIncome::calculateZCBPV(r0, sigma, specialK, rbar, cashflow[i], time[i], steps);
    }

    cout << "The bond price at time 0 is: "<< pv << endl;
}

void RunQn1c(double r0, double sigma, double specialK, double rbar){
    // use explict formula for the underlying bond price
    double FV = 1000;
    double strike = 980;
    double T = 3.0/12;
    double S = 0.5;
    int tradingDays = 252;
    int simNum = 5000;
    int steps = int(floor(T * tradingDays));

    MatrixXd rMat(simNum, steps);
    rMat = MatrixXd::Zero(simNum, steps);

    VectorXd bigR(simNum);
    bigR = VectorXd::Zero(simNum);
    FixedIncome::getRPathAndBigRVasicek(r0, sigma, specialK, rbar, T, steps, rMat, bigR, simNum);


    // the idea is to calcuate the P(T, S) the price at T, which is the exercise price, and S which is the maturity price
    // r_t is the interset rate factor in the explicit formula, it has to be the interest rate at T
    double optionPrice = 0;
    for(int i =0; i< simNum; i++){
        double bondPrice = FixedIncome::getBondPrice(FV, rMat(i, steps - 1), specialK, sigma, rbar, S, T);
        double payoff = Mutils::max(bondPrice - strike, 0); // payoff = bondprice with explicit formula - strike
        optionPrice += payoff/exp(bigR(i));
    }

    cout << "European Call on the bond is:  "<< optionPrice/simNum << endl;
}

void RunQn1d(double r0, double sigma, double specialK, double rbar){
    double strike = 980;
    double T = 0.25;
    double S = 4;
    int tradingDays = 252;
    int simNum = 100;

    int steps = int(floor(T * tradingDays));
    int M = 100;

    vector<double> cashflow{ 30, 30, 30, 30, 30, 30, 30, 1030};
    vector<double> time{0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75};

    // initialize matrix to calculate the Interest Rate at time T
    MatrixXd rMat(simNum, steps);
    rMat = MatrixXd::Zero(simNum, steps);

    VectorXd bigR(simNum);
    bigR = VectorXd::Zero(simNum);
    FixedIncome::getRPathAndBigRVasicek(r0, sigma, specialK, rbar, T, steps, rMat, bigR, simNum);

    // At interest rate at T = 0.25
    VectorXd rAtT(simNum);
    rAtT = rMat.col(steps - 1);

    // The outer for loop generates a vector starting at T with interest rate rAtT(i)
    // Then generate the path from T to S
    double priceAt0 = 0;
    for(int i =0; i< simNum; i++){

        // This is the function to calculate the bond price at the specific interest at time T
        int numCF =8;
        // treat each dividend payment as a zero coupon bond,
        double bondPVatT = 0;
        for(int j=0; j< numCF; j++){
            bondPVatT += FixedIncome::calculateZCBPV(rAtT(i), sigma, specialK, rbar, cashflow[j], time[j], steps);
        }
        //##############################################################

        double p_i_T_S =0;

        for(int j=0; j< M; j++){
            p_i_T_S += Mutils::max(bondPVatT - strike, 0)/exp(-rAtT(i) * T);
        }

        // bond price at T with maturity at S;
        p_i_T_S /= M;

        priceAt0 += p_i_T_S /exp(bigR(i));
    }
    cout << "European Call on the bond is: "<< priceAt0/simNum<<endl;

}

void RunQn2a(double r0, double sigma, double specialK, double rbar){
    double FV = 1000;
    double strike = 980;
    double T = 0.5;
    double S = 1.0;
    double t =0;

    int tradingDays = 252;
    int simNum = 1000;
    int steps = int(floor(T * tradingDays));

    MatrixXd rMat(simNum, steps);
    rMat = MatrixXd::Zero(simNum, steps);

    VectorXd bigR(simNum);
    bigR = VectorXd::Zero(simNum);

    // get the Interest rate process from 0 to 0,5
    FixedIncome::getRPathAndBigRCIR(r0, sigma, specialK, rbar, T, steps, rMat, bigR, simNum);

    double A_t_T = 0.0;
    double B_t_T = 0.0;
    FixedIncome::getfunctionAandB(specialK, sigma, rbar, t, T, A_t_T, B_t_T);


    double optionPrice = 0;
    for(int i =0; i< simNum; i++){
        double bondPrice = FV * A_t_T * exp(-B_t_T * rMat(i, steps - 1));
        double payoff = Mutils::max(bondPrice - strike, 0); // payoff = bondprice with explicit formula - strike
        optionPrice += payoff/exp(bigR(i));
    }

    cout << "European Call on the bond is:  "<< optionPrice/simNum << endl;

}


//
//void RunQn2b(double r0, double sigma, double specialK, double rbar){
//    double FV = 1000;
//    double strike = 980;
//    double T = 0.5;
//    double S = 1.0;
//    double tau = S - T;
//    double t = 0;
//
//    double A_T_S = 0.0;
//    double B_T_S = 0.0;
//    FixedIncome::getfunctionAandB(specialK, sigma, rbar, tau, A_T_S, B_T_S);
//
//    double theta = sqrt(specialK * specialK + 2 * sigma * sigma);
//    double phi = (2 * theta)/(sigma * sigma * (exp(theta * (T-t)) - 1));
//    double psi = (specialK + theta)/(sigma * sigma);
//    double r_star = log(A_T_S/strike)/B_T_S;
//
//    // get P(t, S)
//    double A_t_S = 0.0;
//    double B_t_S = 0.0;
//    FixedIncome::getfunctionAandB(specialK, sigma, rbar, S, A_t_S, B_t_S);
//    double p_t_S = FV * A_t_S * exp(-B_t_S * r0);
//
//    // get P(t, T)
//    double A_t_T_1 = 0.0;
//    double B_t_T_1 = 0.0;
//    FixedIncome::getfunctionAandB(specialK, sigma, rbar, T, A_t_T_1, B_t_T_1);
//    double p_t_T = FV * A_t_T_1 * exp(-B_t_T_1 * r0);
//
//    double chisq_1_x = 2 * r_star * (phi + psi + B_t_S);
//    double chisq_1_p = (4 * specialK * rbar)/(sigma * sigma);
//    double chisq_1_q = (2 * phi * phi * r0 * exp(theta * T))/ (phi + psi + B_t_S);
//
//    double chisq_2_x = 2 * r_star * (phi + psi);
//    double chisq_2_p = (4 * specialK * rbar)/(sigma * sigma);
//    double chisq_2_q = (2 * phi * phi * r0 * exp(theta * T)) /(phi + psi);
//
//
//    double chisq_1_cdf;
//    double chisq_2_cdf;
//
//    double call_option = p_t_S * chisq_1_cdf - strike * p_t_T * chisq_2_cdf;
//
//}
//



int main() {
    double r0 = 0.05;
    double sigma = 0.18;
    double specialK = 0.82;
    double rbar = 0.05;

    cout << "######################################## Qn1 a ##################################################" <<endl;
    //RunQn1a(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn1 b ##################################################" <<endl;
    //RunQn1b(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn1 c ##################################################" <<endl;
    //RunQn1c(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn1 d ##################################################" <<endl;
    //RunQn1d(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;


    r0 = 0.05;
    sigma = 0.12;
    specialK = 0.92;
    rbar = 0.055;

    cout << "######################################## Qn2 a ##################################################" <<endl;
    RunQn2a(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn2 a ##################################################" <<endl;
    //RunQn2b(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

//    cout <<stats::pnorm(1,0,1)<<endl;



    return 0;
}