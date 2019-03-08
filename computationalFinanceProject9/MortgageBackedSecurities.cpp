//
// Created by paragonhao on 6/3/19.
//

#include <cmath>
#include "MortgageBackedSecurities.h"
#include "Mutils.h"
#include <random>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

double MortgageBackedSecurities::getBurnOutRate(const double & pv0, const double & pvt_minus_1){
    return 0.3 + 0.7 * (pvt_minus_1/pv0);
}

double MortgageBackedSecurities::getInterestRate(const double &annualizedMortgageRate,
                                                 const double &tenYearRateT_minus_1){
    return 0.28 + 0.14 * atan(-8.57 + 430 * (annualizedMortgageRate - tenYearRateT_minus_1));
}

// t is the month
double MortgageBackedSecurities::getSeasoningSG(const double &t){
    return Mutils::min(1, t/30.0);
}

double MortgageBackedSecurities::getSeasonalitySY(const int &month_index){
    double val = 0;
    switch (month_index)
    {
        case 0: val = 0.94;
            break;
        case 1: val = 0.76;
            break;
        case 2: val = 0.74;
            break;
        case 3: val = 0.76;
            break;
        case 4: val = 0.95;
            break;
        case 5: val = 0.98;
            break;
        case 6: val = 0.92;
            break;
        case 7: val = 1.10;
            break;
        case 8: val = 1.18;
            break;
        case 9: val = 1.22;
            break;
        case 10: val = 1.23;
            break;
        case 11: val = 0.98;
            break;
        default: printf("wrong input");
            break;
    }
    return val;
}


/* Get CPR at time t
 *
 * @param:
 * double pv0 : PV of loan at time 0
 * double pvt_minus_1: previous pv at t -1
 * double t: from 1 to 360, for each month
 * double R: annualize rate R (wac), constant
 * double rt: interest rate at time t
 * int window: the 10 year interest rate window, set to 120 (month)
 * double T: number of years, 30 year in this case
 * int steps: number of steps, 361 in total in this case,
 * double specialK: ki, constant
 * double rbar: rbar, constant
 * double sigma: sigma, constant,
 * int monthindex: 1 to 12;
 * */
double MortgageBackedSecurities::getCPR_t(const double &pv0, const double &pvt_minus_1, const double &t,
                                          const double &R, const double &rt, const int &window, const double &T,
                                          const int &steps, const double &specialK,
                                          const double &rbar, const double &sigma, const int & monthIndex){
    double curr_r = rt;
    double delta_t = 10 / steps;

    double r_t_minus_1 = 0;
    // rt is the starting pos for the interest rate process
    srand(time(0));
    int seed = rand();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);

    for(int i=0; i<window; i++){
        double z = distribution(generator);
        double r_next = curr_r + specialK * (rbar - curr_r) * delta_t + sigma * sqrt(abs(curr_r)) * sqrt(delta_t) * z;
        r_t_minus_1 += r_next ;
        curr_r = r_next;
    }
    r_t_minus_1 /= window;

    double CPR_t = getInterestRate(R, r_t_minus_1) * getBurnOutRate(pv0, pvt_minus_1) * getSeasoningSG(t) * getSeasonalitySY(monthIndex);

    return CPR_t;

}

/* Get interest rate process over T years
 *
 * @param:
 * r0 : starting interest rate
 * double sigma: sigma
 * double specialK: ki
 * double rbar: r bar
 * double T: Period during which the process exists
 * int steps: number of steps
 * MatrixXd &rMat: initial stock price of x
 * int simNum: number of simulations
 * */
void MortgageBackedSecurities::getRPathCIRModel(double r0, double sigma, double specialK, double rbar, double T, int steps,
                                        MatrixXd &rMat, int simNum){

    double delta_t = T/steps;

    srand(time(0));
    int seed = rand();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);


    // clear Matrix rMat and bigR
    rMat = MatrixXd::Zero(simNum, steps);

    // initialize r[0]  to be r0
    for(int i =0; i < simNum ;i++){
        rMat(i,0) = r0;
    }

    // Discretization of the r matrix
    for(int i =0; i< simNum; i++){
        for(int j =1; j< steps; j++){
            double z = distribution(generator);
            rMat(i, j) =  rMat(i, j-1) + specialK * (rbar -  rMat(i, j-1)) * delta_t + sigma * sqrt(abs(rMat(i, j-1))) * sqrt(delta_t) * z;
        }
    }
}