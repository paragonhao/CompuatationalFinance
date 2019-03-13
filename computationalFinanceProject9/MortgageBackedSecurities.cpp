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
                                          const double &R, const double &rt, const int & monthIndex){

    double CPR_t = getInterestRate(R, rt) * getBurnOutRate(pv0, pvt_minus_1) * getSeasoningSG(t) * getSeasonalitySY(monthIndex);

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

//    srand(time(0));
    int seed = 12345678;
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0.0,1.0);

    // clear Matrix rMat and bigR
    rMat = MatrixXd::Zero(simNum, steps);

    // Discretization of the r matrix
    for(int i =0; i< simNum; i++){
        rMat(i,0) = r0;
        for(int j =1; j< steps; j++){
            double z = distribution(generator);
            rMat(i, j) =  rMat(i, j-1) + specialK * (rbar -  rMat(i, j-1)) * delta_t + sigma * sqrt(abs(rMat(i, j-1))) * sqrt(delta_t) * z;
        }
    }
}



double MortgageBackedSecurities::getNumerixPrepaymentModel(double WAC, double pv0, double r0,
                                                         double specialK, double rbar, double sigma,
                                                         double time, double OAS, bool isOAS){
    cout << "Calculating MBS"<<endl;
    double x = isOAS?OAS:0;
    cout <<"Trying OAS Spread: " <<x<<endl;

    int interval = 12; // 12 month a year or change to 360 days per year
    int window = 120; // 10 year = 120 months to calculate 10 year rate
    double r = WAC/12;

// step 1: generating interest rate process for 30 years, this r path would be used as cash flow discount factor
    int simNum = 5000;
    auto steps = int(time * interval + 1); // taking each month  should be 361, including the initial r0;
    double delta_t = time/steps;

    MatrixXd rMat;
    rMat = MatrixXd::Zero(simNum, steps);
    MortgageBackedSecurities::getRPathCIRModel(r0, sigma, specialK, rbar, (time + 10), (steps + 10 * interval), rMat, simNum);

    MatrixXd r10Mat;
    r10Mat = MatrixXd::Zero(simNum, steps);
    for(int i=0; i< simNum;i++){
        for(int j=0; j<steps; j++){
            for(int l=j; l< window + j; l++){
                r10Mat(i,j) += rMat(i,l) * delta_t;
            }
//            r10Mat(i,j) /= window;
        }
    }

    double mbsPrice = 0;

    for(int sim =0; sim< simNum;sim++){

        double pv_t_minus_1 = pv0;
        double pv_this_iter = 0;

        for(int i=1; i<steps; i++){


            double CPR_t = MortgageBackedSecurities::getCPR_t(pv0, pv_t_minus_1, i, WAC, r10Mat(sim, i), (i+1)%12);
            double temp_1 = (1.0/ (1.0 - pow((1.0 + r),(i - 1.0 - steps)))) - 1.0;
            double temp_2 = pow((1 - CPR_t), 1.0/12);
            double tpp_t = pv_t_minus_1 * r * temp_1 + (pv_t_minus_1 - pv_t_minus_1 * r * temp_1) * (1 - temp_2);
            double pv_t = pv_t_minus_1 - tpp_t;
            double c_t = tpp_t + pv_t * r; //in total there should be 360 cashflows

            double R =0;
            for(int m=1; m<=i; m++){
                R += (rMat(sim,m) + x) * delta_t;
            }
            pv_this_iter += exp(-R) * c_t;

            pv_t_minus_1 = pv_t;
        }
        mbsPrice += pv_this_iter;
    }

    cout <<"rbar = "<<rbar <<", sigma = "<< sigma <<", k = " << specialK <<", Price = "<< mbsPrice/double(simNum)<< endl;
    return mbsPrice/double(simNum);
}

