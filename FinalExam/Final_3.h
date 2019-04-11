//
// Created by paragonhao on 20/3/19.
//

#ifndef FINALEXAM_FINAL_3_H
#define FINALEXAM_FINAL_3_H


#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <iostream>
#include "Mutils.h"

using namespace Eigen;
using namespace std;

class Qn3 {

    public:


    static int RunQ3() {

        cout << "Running Qn 3"<<endl;
        double k = 100;
        double T = 5.0;
        double s0 = 100;
        double r = 0.05;
        double sigma = 0.35;
        int nSteps = 101;
        int nSim = 10000;
        double delta_t = T/100.0;

        // first, generate stock price process
        int seed = 12345678;
        default_random_engine generator(seed);
        normal_distribution<double> distribution(0.0,1.0);

        MatrixXd priceProcess = MatrixXd::Zero(nSim, nSteps);

        int half = nSim/2;
        for(int i =0; i< half; i++){
            priceProcess(i,0) = s0;
            priceProcess(half +i, 0) =s0;
            for(int j =1; j< nSteps; j++){
                double z = distribution(generator);
                priceProcess(i, j) =  priceProcess(i, j-1) * exp((r - 0.5 * sigma * sigma) * delta_t + sigma * z * sqrt(delta_t));
                priceProcess(half + i, j) =  priceProcess(half + i, j-1) * exp((r - 0.5 * sigma * sigma) * delta_t + sigma * (-z) * sqrt(delta_t));
            }
        }


        // second generate upper and lower bound
        VectorXd myTime = VectorXd::Zero(nSteps);

        for(int i=1;i< nSteps;i++){
            myTime(i) = delta_t + myTime(i-1);
        }

        // create upper lower vector
        VectorXd upper = VectorXd::Zero(nSteps);
        VectorXd lower = VectorXd::Zero(nSteps);

        for(int i=0; i< nSteps; i++){
            lower(i) = 50 * exp(0.138629 * myTime(i));
            upper(i) = 200 - 50 * exp(0.138629 * myTime(i));
        }

        // find the payoff
        int counterForLower =0;

        double payOff = 0;

        for(int i=0;i< nSim;i++){

            for(int j=0; j< nSteps;j++){
                double currPrice = priceProcess(i,j); //current price

                if(currPrice > upper(j)){ // current price > upper bound, ex like call payoff
                    payOff +=Mutils::max(currPrice - k, 0);
                    break; // once pay off is collected, break out of this path
                }else if(currPrice < lower(j)){  // current price < lower bound, ex like call payoff
                    payOff += Mutils::max(k - currPrice, 0);
                    counterForLower ++ ; // count the number of times price hit lower bound first
                    break;
                }
            }
        }
        cout <<"pay off: " << (payOff/nSim) * exp(-r * T)<<endl;
        cout << "conditional probability: "<< (counterForLower*1.0)/nSim<<endl;

        return 0;
    }

};
#endif //FINALEXAM_FINAL_3_H
