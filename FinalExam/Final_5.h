//
// Created by paragonhao on 20/3/19.
//

#ifndef FINALEXAM_FINAL_5_H
#define FINALEXAM_FINAL_5_H


#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <iostream>
#include "Mutils.h"

using namespace Eigen;
using namespace std;

class Qn5 {

    public:

    static void RunQ5() {
        cout << "Running Qn 5"<<endl;
        double s0 = 6000;
        double rho = -0.25;
        double r = 0.05;
        double q = 0;
        double rf = 0.04;
        double sigma1 = 0.1;
        double sigma2 = 0.15;
        double gamma = -0.04;
        double lamda = 1.5;
        double k = 60;
        double T = 1.0;
        double E0 = 0.0096;
        int nSim = 10000;
        int nSteps = 100;
        double delta_t = T/nSteps;

        MatrixXd dW;
        Mutils::generateWProcessMat(nSim, nSteps, delta_t, dW, 123456);

        MatrixXd dW1;
        Mutils::generateWProcessMat(nSim, nSteps, delta_t, dW1, 654321);

        MatrixXd dB;

        dB = MatrixXd::Zero(nSim, nSteps);

        for(int i=0; i< nSim; i++){
            for(int j=0; j<nSteps; j++){
                dB(i,j) = rho * dW(i,j) + sqrt(1 - rho * rho) * dW1(i,j);
            }
        }

        int currSteps = nSteps + 1;
        MatrixXd S_t = MatrixXd::Zero(nSim, currSteps);
        MatrixXd E_t = MatrixXd::Zero(nSim, currSteps);

        for(int i=0; i< nSim; i++){
            S_t(i, 0) = s0;
            E_t(i, 0) = E0;
        }

        std::default_random_engine rpois(12345678);
        std::poisson_distribution<int> distribution(lamda * delta_t);

        // Use dW and dB
        for(int i=0; i< nSim; i++){
            for(int j=0; j<nSteps; j++){
                double pois = distribution(rpois);
                S_t(i, j+1) = S_t(i, j) + S_t(i, j) * ((r - q) * delta_t + sigma1 * dW(i, j) + gamma * pois);
                E_t(i, j+1) = E_t(i, j) + E_t(i, j) * delta_t * (r - rf) + sigma2 * E_t(i, j) * dB(i, j);
            }
        }

        double payOff =0;
        for(int i=0; i< nSim; i++){
            payOff += Mutils::max(S_t(i, currSteps-1) * E_t(i, currSteps-1) - k, 0);
        }

        cout << "Pay off is: " <<(payOff/nSim) * exp(-r*T)<<endl;
    }

};
#endif //FINALEXAM_FINAL_5_H
