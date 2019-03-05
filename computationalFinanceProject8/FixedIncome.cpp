//
// Created by paragonhao on 3/3/19.
//

#include "FixedIncome.h"
#include "RandomGenerator.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <ctime>

using namespace Eigen;
using namespace std;

double FixedIncome::calculateZCBPV(double r0, double sigma, double specialK, double rbar, double FV, double T, int steps){

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
            rMat(i, j) =  rMat(i, j-1) + specialK * (rbar -  rMat(i, j-1)) * delta_t + sigma * sqrt(delta_t) * stdNormArray[stdNormCounter++];
            bigR(i) += rMat(i, j) * delta_t;
        }
    }

    double priceAt0 = 0;
    for(int i =0; i< simNum;i++){
        priceAt0 += FV/exp(bigR(i));
    }
    return priceAt0/simNum;
}

double FixedIncome::getInterestRateSim(double r0, double sigma, double specialK, double rbar, double T, int steps, double payoff) {

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

    double optionPrice = 0;
    for(int i =0; i< simNum; i++){
        optionPrice += payoff/exp(bigR(i));
    }
    return optionPrice/simNum;
}


double FixedIncome::getBondPrice(double FV, double rt, double specialK, double sigma, double rbar, double T, double t){

    double B = (1/specialK) * (1 - exp(-specialK * (T - t)));

    double A = exp((rbar - (sigma * sigma) / ( 2 * specialK * specialK)) * (B - (T - t)) - (sigma * sigma * B * B)/(4 * specialK));

    double bondPrice = A * exp(-B * rt) * FV;

    return bondPrice;
}

void FixedIncome::getRPathAndBigR(double r0, double sigma, double specialK, double rbar, double T, int steps, MatrixXd & rMat, VectorXd & bigR){

    double delta_t = T/steps;
    int simNum = 1000;
    // create a 126 * 1000 array of std normal distribution, total should be 126000 number of Zs;
    int N = (steps - 1) * simNum;
    auto * stdNormArray = new double[N];

    srand(std::time(0)); //use current time as seed for random generator
    int seed = rand();

    stdNormArray = RandomGenerator::boxmuller(RandomGenerator::runif(N, seed), N);

    // clear Matrix rMat and bigR
    rMat = MatrixXd::Zero(simNum, steps);

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
}