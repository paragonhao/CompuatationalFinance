//
// Created by paragonhao on 3/3/19.
//

#include "FixedIncome.h"
#include "RandomGenerator.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <ctime>
#include <random>
#include <iostream>

using namespace Eigen;
using namespace std;

double FixedIncome::calculateZCBPV(double r0, double sigma, double specialK, double rbar, double FV, double T, int steps){

    int simNum = 100; // rows

    MatrixXd rMat(simNum, steps);
    rMat = MatrixXd::Zero(simNum, steps);

    VectorXd bigR(simNum);
    bigR = VectorXd::Zero(simNum);

    FixedIncome::getRPathAndBigRVasicek(r0, sigma, specialK, rbar, T, steps, rMat, bigR, simNum);

    double priceAt0 = 0;
    for(int i =0; i< simNum;i++){
        priceAt0 += FV/exp(bigR(i));
    }
    return priceAt0/simNum;
}

double FixedIncome::getBondPrice(double FV, double rt, double specialK, double sigma, double rbar, double T, double t){

    double B = (1/specialK) * (1 - exp(-specialK * (T - t)));

    double A = exp((rbar - (sigma * sigma) / ( 2 * specialK * specialK)) * (B - (T - t)) - (sigma * sigma * B * B)/(4 * specialK));

    double bondPrice = A * exp(-B * rt) * FV;

    return bondPrice;
}

void FixedIncome::getRPathAndBigRVasicek(double r0, double sigma, double specialK, double rbar, double T, int steps,
                                         MatrixXd &rMat, VectorXd &bigR, int simNum){

    double delta_t = T/steps;

    int seed = rand();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);


    // clear Matrix rMat and bigR
    rMat = MatrixXd::Zero(simNum, steps);

    bigR = VectorXd::Zero(simNum);

    // initialize r[0]  to be r0
    for(int i =0; i < simNum ;i++){
        rMat(i,0) = r0;
    }

    // Discretization of the r matrix
    for(int i =0; i< simNum; i++){
        for(int j =1; j< steps; j++){
            double z = distribution(generator);
            rMat(i, j) =  rMat(i, j-1) + specialK * (rbar -  rMat(i, j-1)) * delta_t + sigma * sqrt(delta_t) * z;
            bigR(i) += rMat(i, j) * delta_t;
        }
    }
}


void FixedIncome::getRPathAndBigRCIR(double r0, double sigma, double specialK, double rbar, double T, int steps,
                                         MatrixXd &rMat, VectorXd &bigR, int simNum){

    double delta_t = T/steps;

    int seed = rand();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);


    // clear Matrix rMat and bigR
    rMat = MatrixXd::Zero(simNum, steps);

    bigR = VectorXd::Zero(simNum);

    // initialize r[0]  to be r0
    for(int i =0; i < simNum ;i++){
        rMat(i,0) = r0;
    }

    // Discretization of the r matrix
    for(int i =0; i< simNum; i++){
        for(int j =1; j< steps; j++){
            double z = distribution(generator);
            rMat(i, j) =  rMat(i, j-1) + specialK * (rbar -  rMat(i, j-1)) * delta_t + sigma * sqrt(rMat(i, j-1)) * sqrt(delta_t) * z;
            bigR(i) += rMat(i, j) * delta_t;
        }
    }
}

void FixedIncome::getfunctionAandB(double specialK, double sigma, double rbar, double t, double T, double &A_t_T, double &B_t_T){
    double h1 = sqrt(specialK * specialK  + 2 * sigma * sigma);
    double h2 = (specialK + h1) * 0.5;
    double h3 = (2 * specialK * rbar) / (sigma * sigma);
    double exp_h1_tau = exp(h1 * (T - t));
    double exp_h2_tau = exp(h2 * (T - t));

    A_t_T = pow((h1 * exp_h2_tau )/((h2 * (exp_h1_tau - 1)) + h1), h3);
    B_t_T = (exp_h1_tau - 1)/((h2 * (exp_h1_tau - 1)) + h1);
}