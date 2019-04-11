//
// Created by paragonhao on 20/3/19.
//

#ifndef FINALEXAM_FINAL_2_H
#define FINALEXAM_FINAL_2_H


#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <iostream>
#include "Mutils.h"

using namespace Eigen;
using namespace std;

class Qn2 {

    public:

    static void HestonModelSimulation(double v_0, double alpha, double beta, double gamma, double s0, double r, double rho, double T, double nSteps, double nSim){
        // genereate a wp process
        double delta_t = T / nSteps;
        MatrixXd dW;
        //    MatrixXd dW;
        Mutils::generateWProcessMat(nSim, nSteps, delta_t, dW, 123456);

        // genereate another wp process
        MatrixXd dW1;
        Mutils::generateWProcessMat(nSim, nSteps, delta_t, dW1, 654321);

        MatrixXd dB;
        dB = MatrixXd::Zero(nSim, nSteps);

        // generate correlated process
        for(int i=0; i<nSim; i++){
            for(int j=0;j<nSteps;j++){
                dB(i,j) = rho * dW(i,j) + sqrt(1 - rho * rho) * dW1(i,j);
            }
        }

        //use dw and dB
        int currSteps = nSteps + 1;
        // generate Stock process and volatility
        MatrixXd S_t = MatrixXd::Zero(nSim, currSteps);
        MatrixXd V_t = MatrixXd::Zero(nSim, currSteps);

        for(int i=0; i<nSim; i++){
            S_t(i, 0) = s0;
            V_t(i, 0) = v_0;
            for(int j=0; j< nSteps; j++){
                S_t(i, j+1) = S_t(i, j) + S_t(i, j) * r * delta_t + sqrt( Mutils::max(V_t(i,j),0)) * S_t(i, j) * dW(i,j);
                V_t(i, j+1) = V_t(i, j) + (alpha + beta * Mutils::max(V_t(i,j),0)) * delta_t + gamma * sqrt(Mutils::max(V_t(i,j),0)) * dB(i,j);
            }
        }

        // calculate the payoff of the option
        double payOff = 0;

        for(int i =0;i<nSim; i++){
            double avg = S_t.row(i).mean();
            payOff += Mutils::max(S_t(i, currSteps - 1) - avg,0);
        }

        cout << (payOff /nSim) * exp(-r * T)<<endl;

    }


    static int RunQ2() {
        cout << "Running Qn 2"<<endl;
        double v_0 = 0.06;
        double alpha = 0.45;
        double beta = -5.105;
        double gamma = 0.25;
        double s0 = 20;
        double r = 0.05;
        double rho = -0.75;
        double T = 2.0;

        int nSteps = 100;
        int nSim = 10000;

        HestonModelSimulation( v_0,  alpha,  beta,  gamma,  s0,  r,  rho,  T,  nSteps,  nSim);
        rho = 0;
        HestonModelSimulation( v_0,  alpha,  beta,  gamma,  s0,  r,  rho,  T,  nSteps,  nSim);
        rho = 0.75;
        HestonModelSimulation( v_0,  alpha,  beta,  gamma,  s0,  r,  rho,  T,  nSteps,  nSim);
        return 0;
    }

};
#endif //FINALEXAM_FINAL_2_H
