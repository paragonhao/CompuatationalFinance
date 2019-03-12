//
// Created by paragonhao on 3/3/19.
//

#include "FixedIncome.h"
#include "RandomGenerator.h"
#include "Mutils.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <ctime>
#include <random>
#include <iostream>

using namespace Eigen;
using namespace std;

double FixedIncome::calculateZCBPV(double r0, double sigma, double specialK, double rbar, double FV, double T, int steps){

    int simNum = 1000; // rows

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

    srand(time(0));
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


void FixedIncome::getExplicitOptionCallPriceCIR(double r0, double sigma, double specialK, double rbar,
                                                double strike, double FV, double T, double S, double t){
    double h1 = sqrt(specialK*specialK + 2 * sigma*sigma);
    double h2 = (specialK + h1) / 2;
    double h3 = 2 * specialK*rbar / (sigma*sigma);
    double rS = r0 * exp(-specialK * (T - t)) + rbar * (1 - exp(-specialK * (T - t)));

    double A_T_S = pow((h1*exp(h2*(S - T))) / (h2*(exp(h1*(S - T)) - 1) + h1), h3);
    double B_T_S = (exp(h1*(S - T)) - 1) / (h2*(exp(h1*(S - T)) - 1) + h1);
    double P_T_S = A_T_S * exp(-B_T_S * rS);

    double A_t_S = pow((h1*exp(h2*(S - t))) / (h2*(exp(h1*(S - t)) - 1) + h1), h3);
    double B_t_S = (exp(h1*(S - t)) - 1) / (h2*(exp(h1*(S - t)) - 1) + h1);
    double P_t_S = A_t_S * exp(-B_t_S * r0);

    double A_t_T = pow((h1*exp(h2*(T - t))) / (h2*(exp(h1*(T - t)) - 1) + h1), h3);
    double B_t_T = (exp(h1*(T - t)) - 1) / (h2*(exp(h1*(T - t)) - 1) + h1);
    double P_t_T = A_t_T * exp(-B_t_T * r0);

    double theta = h1;
    double phi = 2 * theta / (sigma*sigma*(exp(theta*(T - t)) - 1));
    double psi = (specialK + theta) / (sigma*sigma);
    double rStar = log(A_T_S /(strike/FV))/ B_T_S;

    double x1 = 2 * rStar * (phi + psi + B_T_S);
    double p1 = 4 * specialK * rbar/(sigma*sigma);
    double q1 = (2 * phi*phi*r0 * exp(theta*(T - t))) / (phi + psi + B_T_S);
    double x2 = 2 * rStar*(phi + psi);
    double p2 = 4 * specialK*rbar / (sigma*sigma);
    double q2 = (2 * phi*phi*r0 * exp(theta*(T - t))) / (phi + psi);


    //these two values were calculated by matlab
    double chisquared1 = 0.2893;
    double chisquared2 = 0.2864;
    double call = FV * P_t_S * chisquared1 - strike * P_t_T * chisquared2;
    cout <<"Call option price is: "<< call << endl;

}

void FixedIncome::getRPathAndBigRCIR(double r0, double sigma, double specialK, double rbar, double T, int steps,
                                         MatrixXd &rMat, VectorXd &bigR, int simNum){

    double delta_t = T/steps;

    srand(time(0));
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

void FixedIncome::getRPathAndBigRGPlusPlusModel(int simNum, double T, int steps, double a, double b, double sigma, double phi, double eta,
                                                double r0, double rho, MatrixXd &rMat, VectorXd & bigR, MatrixXd &xMat, MatrixXd &yMat){

    double delta = T/steps;

    srand(time(0));
    int seed = rand();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);

    for(int i =0; i< simNum; i++){
        for(int j =1; j< steps; j++){
            double z1 = distribution(generator);
            double z2 = rho * z1 + sqrt(1-(rho * rho)) * distribution(generator);
            xMat(i, j) = xMat(i, j-1) + (-a * xMat(i, j-1)) * delta + sigma * z1 * sqrt(delta);
            yMat(i, j) = yMat(i, j-1) + (-b * yMat(i, j-1)) * delta + eta * z2 * sqrt(delta);
            rMat(i,j) = xMat(i, j) + yMat(i, j) + phi;
            bigR(i) += rMat(i, j) * delta;

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
