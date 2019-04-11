//
// Created by paragonhao on 20/3/19.
//

#ifndef FINALEXAM_FINAL_4_H
#define FINALEXAM_FINAL_4_H


#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <iostream>
#include "Mutils.h"

using namespace Eigen;
using namespace std;

class Qn4 {

    public:

    static void RunQ4() {

        cout << "Running Qn 4"<<endl;
        double r0 = 0.05;
        double alpha = 0.36;
        double beta = -5.86;
        double sigma = 0.36;
        double gamma = 2.0;
        double k = 9800;
        double T = 0.5; // Option expiration
        double S = 1.0; // bond maturity
        int nSteps = 100; //cols
        int nSims = 10000; //rows
        double delta_t = S/nSteps; // delta is 0.01
        double FV = 10000;


        // first generate the interest rate process
        // generate the WP process
//        MatrixXd dW;
//        dW =  MatrixXd::Zero(nSims, nSteps);

        int seed = 654321;
        default_random_engine generator(seed);
        normal_distribution<double> distribution(0.0, 1.0);

//        Mutils::generateWProcessMat(nSims, nSteps, delta_t, dW, seed);

        MatrixXd r_t = MatrixXd::Zero(nSims, nSteps);

        for(int i =0; i< nSims; i++){
            r_t(i,0)= r0;
            for(int j=1; j< nSteps; j++){
                double z = distribution(generator);
                r_t(i,j) = r_t(i, j-1) + (alpha + beta * r_t(i, j-1)) * delta_t + sigma * pow(r_t(i, j-1), gamma) * z * delta_t;
            }
        }

        // second, find the interest rate at time T and S,
        VectorXd r_zero_To_T = VectorXd::Zero(nSims);
        VectorXd r_T_To_S = VectorXd::Zero(nSims);

        int zero_To_T_num = nSteps/2;
        int T_To_S_num = nSteps;

        // for each simulation path, use delta * interest rate
        for(int i =0; i< nSims; i++){
            double curr_r = 0;
            for(int j =0; j< zero_To_T_num; j++){

                curr_r += r_t(i,j) * delta_t;
            }
            r_zero_To_T(i) = curr_r;
        }

        //find out interest from T to S
        for(int i =0;i<nSims; i++){
            double curr_r =0;
            for(int j=zero_To_T_num; j< T_To_S_num;j++){
                curr_r += r_t(i,j) * delta_t;
            }
            r_T_To_S(i) = curr_r;
        }


        // next, calculate the price of the bond at time T,using r_T_To_S to discount
        VectorXd p_t_s = VectorXd::Zero(nSims);

        for(int i =0; i< nSims; i++){
            p_t_s(i) = FV * exp(-r_T_To_S(i));
        }

        double payOff_V = 0;
        for(int i=0; i< nSims; i++){
            payOff_V += Mutils::max(k - p_t_s(i),0) * exp(-r_zero_To_T(i));
        }

        cout << payOff_V/nSims<<endl;
    }

};
#endif //FINALEXAM_FINAL_4_H
