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
    int simNum = 1000;

    int steps = int(floor(T * tradingDays));
    int M = 1000;

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
            p_i_T_S += Mutils::max(bondPVatT - strike, 0);
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

    int tradingDays = 360;
    int simNum = 5000;
    int steps = int(floor(T * tradingDays));

    MatrixXd rMat(simNum, steps);
    rMat = MatrixXd::Zero(simNum, steps);

    VectorXd bigR(simNum);
    bigR = VectorXd::Zero(simNum);

    // get the Interest rate process from 0 to 0,5
    FixedIncome::getRPathAndBigRCIR(r0, sigma, specialK, rbar, T, steps, rMat, bigR, simNum);

    double A_t_T = 0.0;
    double B_t_T = 0.0;
    FixedIncome::getfunctionAandB(specialK, sigma, rbar, T, S, A_t_T, B_t_T);

    double optionPrice = 0;
    for(int i =0; i< simNum; i++){
        double bondPrice = FV * A_t_T * exp(-B_t_T * rMat(i, steps - 1));
        double payoff = Mutils::max(bondPrice - strike, 0); // payoff = bondprice with explicit formula - strike

        optionPrice += payoff/exp(bigR(i));
    }

    cout << "European Call on the bond is:  "<< optionPrice/simNum << endl;

}

void RunQn2b(double r0, double sigma, double specialK, double rbar){
    double FV = 1000;
    double strike = 980;
    double t = 0;
    double T = 0.5;
    double S = 1;
    FixedIncome::getExplicitOptionCallPriceCIR(r0,  sigma,  specialK,  rbar,  strike, FV, T,  S,  t);

}


void RunQn3(double r0, double sigma){
    // expiry time and maturity time
    double T = 0.5;
    double S = 1;
    double strike = 985;
    int steps = int(floor(T * 360)) + 1;
    int simNum =1000;
    double FV =1000;

    // variables
    double rho = 0.7;
    double a = 0.1;
    double b = 0.3;
    double eta = 0.08;
    double phi = 0.03;

    MatrixXd rMat(simNum, steps);
    rMat = MatrixXd::Zero(simNum, steps);

    for(int i = 0; i< simNum; i++){
        rMat(i, 0) = r0;
    }

    MatrixXd xMat(simNum, steps);
    xMat = MatrixXd::Zero(simNum, steps);

    MatrixXd yMat(simNum, steps);
    yMat = MatrixXd::Zero(simNum, steps);

    VectorXd bigR(simNum);
    bigR = VectorXd::Zero(simNum);

    FixedIncome::getRPathAndBigRGPlusPlusModel(simNum, T, steps, a, b, sigma, phi,  eta,  r0,  rho, rMat, bigR, xMat, yMat);

    double optionPrice =0;

    for(int i=0; i< simNum; i++) {
        MatrixXd rMatT_S(simNum, steps);
        rMatT_S = MatrixXd::Zero(simNum, steps);

        // using one fixed r as starting point for this iteration
        for (int j = 0; j < simNum; j++) {
            rMatT_S(j, 0) = rMat(i, steps - 1);
        }

        MatrixXd xMatT_S(simNum, steps);
        xMatT_S = MatrixXd::Zero(simNum, steps);
        for (int j = 0; j < simNum; j++) {
            xMatT_S(j, 0) = xMat(i, steps - 1);
        }

        MatrixXd yMatT_S(simNum, steps);
        yMatT_S = MatrixXd::Zero(simNum, steps);
        for (int j = 0; j < simNum; j++) {
            yMatT_S(j, 0) = yMat(i, steps - 1);
        }

        VectorXd bigRT_S(simNum);
        bigRT_S = VectorXd::Zero(simNum);


        FixedIncome::getRPathAndBigRGPlusPlusModel(simNum, (S - T), steps, a, b, sigma, phi, eta, r0, rho, rMat, bigRT_S,
                                                   xMatT_S, yMatT_S);
        double payoff = 0;
        for (int m = 0; m < simNum; m++) {
            double bondPrice = FV * exp(-bigRT_S(m));
            payoff += Mutils::max(strike - bondPrice, 0);
        }

        payoff /= simNum;

        optionPrice += payoff * exp(-bigR(i));

    }
    cout << "Option Price is: "<<optionPrice / simNum << endl;
}



int main() {
    double r0 = 0.05;
    double sigma = 0.18;
    double specialK = 0.82;
    double rbar = 0.05;

    cout << "######################################## Qn1 a ##################################################" <<endl;
    RunQn1a(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn1 b ##################################################" <<endl;
    RunQn1b(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn1 c ##################################################" <<endl;
    RunQn1c(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn1 d ##################################################" <<endl;
    RunQn1d(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;


    r0 = 0.05;
    sigma = 0.18;
    specialK = 0.92;
    rbar = 0.055;

    cout << "######################################## Qn2 a ##################################################" <<endl;
    RunQn2a(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    cout << "######################################## Qn2 b ##################################################" <<endl;
    RunQn2b(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    r0 = 0.03;
    sigma = 0.03;

    cout << "######################################## Qn3 ##################################################" <<endl;
    RunQn3(r0, sigma);
    cout << "###############################################################################################" <<endl;

    return 0;
}