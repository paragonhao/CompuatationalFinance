//
// Created by paragonhao on 20/3/19.
//

#ifndef FINALEXAM_FINAL_1_H
#define FINALEXAM_FINAL_1_H


#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <cmath>
#include <iostream>
#include "Mutils.h"

using namespace Eigen;
using namespace std;

class Qn1 {

    public:
    static int RunQ1() {
        cout << "Running Qn 2"<<endl;
        //generate
        int nSim =10000;

        int seed = 12345678;
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distribution(0.0,1.0);

        double x = 1.1;
        double payoff = 0;

        for(int i =0; i< nSim; i++){
            double S_k = 0;
            int pos = 0;

            while(S_k < x){
                double num = distribution(generator);
                S_k += num;
                pos++;
            }
            //pos is the value of k in this simulation
            payoff += Mutils::max(4.54 - pos, 0);
        }

        cout << payoff/nSim<<endl;

        return 0;
    }

};
#endif //FINALEXAM_FINAL_1_H
