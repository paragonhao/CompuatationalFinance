#include <iostream>
#include "MortgageBackedSecurities.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void RunQn1a(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double time){

    int interval = 12; // 12 month a year or change to 360 days per year
    int window = 120; // 10 year = 120 months to calculate 10 year rate
    double r = r0/12;

// step 1: generating interest rate process for 30 years, this r path would be used as cash flow discount factor
    int simNum = 1000;
    auto steps = int(time * interval + 1); // taking each month  should be 361, including the initial r0;
    MatrixXd rMat;
    rMat = MatrixXd::Zero(simNum, steps);
    MortgageBackedSecurities::getRPathCIRModel(r0, sigma, specialK, rbar, time, steps, rMat, simNum);

    double delta_t = time/steps;

    double mbsPrice = 0;

    for(int sim =0; sim< simNum;sim++){

        double pv_t_minus_1 = pv0;
        double pv_this_iter = 0;

        for(int i=1; i<steps; i++){

            double temp_1 = (1/ (1 - pow((1 + r),(i - 1 - steps)))) - 1;
            double CPR_t = MortgageBackedSecurities::getCPR_t(pv0,pv_t_minus_1, i,WAC,rMat(sim,i), window, time, steps, specialK, rbar, sigma, (i+1)%12);
            double tpp_t = pv_t_minus_1 * r * temp_1 + (pv_t_minus_1 - pv_t_minus_1 * r * temp_1) * (1 - pow(1 - CPR_t, 1/12));
            double pv_t = pv_t_minus_1 - tpp_t;
            double c_t = tpp_t + pv_t * r; //in total there should be 360 cashflows

            double R =0;
            for(int m=1; m<=i; m++){
                R += rMat(sim,m) * delta_t;
            }
            pv_this_iter += exp(-R) * c_t;

            pv_t_minus_1 = pv_t;
        }

        mbsPrice += pv_this_iter;

    }

    cout << mbsPrice/simNum<< endl;
}



int main() {
    double WAC = 0.08;
    double pv0 = 100000;
    double r0 = 0.078;
    double specialK = 0.9;
    double rbar = 0.08;
    double sigma = 0.12;
    double time = 30; //num of years

    RunQn1a( WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time);

    return 0;
}