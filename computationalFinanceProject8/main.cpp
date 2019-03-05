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
    int steps = int(floor(T * tradingDays));

    double bondPrice = FixedIncome::getBondPrice(FV, r0, specialK, sigma, rbar, S, T);
    double payoff = Mutils::max(bondPrice - strike, 0);


    //Generate R from 0 to T(0.25) first
    double delta_t = T/steps;
    int simNum = 1000; // rows

    // create a 126 * 1000 array of std normal distribution, total should be 126000 number of Zs;
    int N = (steps - 1) * simNum;
    auto * stdNormArray = new double[N];

    srand(std::time(0)); //use current time as seed for random generator
    int seed = rand();

    stdNormArray = RandomGenerator::boxmuller(RandomGenerator::runif(N, seed), N);

    MatrixXd rMat(simNum, steps);
    rMat = MatrixXd::Zero(simNum, steps);

    VectorXd bigR(simNum);
    bigR = VectorXd::Zero(simNum);

    // initialize r[0]  to be r0
    for(int i =0; i < simNum ;i++){
        rMat(i,0) = r0;
    }

    // Discretization of the r matrix
    int stdNormCounter = 0;
    for(int i =0; i< simNum; i++){
        for(int j =1; j< steps; j++){
            rMat(i, j) =  abs(rMat(i, j-1)) + specialK * (rbar -  abs(rMat(i, j-1))) * delta_t + sigma * sqrt(delta_t) * stdNormArray[stdNormCounter++];
            bigR(i) += rMat(i, j) * delta_t;
        }
    }
    //#######################################################################################################


    double optionPrice = 0;
    for(int i =0; i< simNum; i++){
        optionPrice += payoff/exp(bigR(i));
    }

    cout << "The option price at time 0 is: "<< optionPrice/simNum << endl;

}


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
    RunQn1c(r0, sigma, specialK, rbar);
    cout << "###############################################################################################" <<endl;

    return 0;
}