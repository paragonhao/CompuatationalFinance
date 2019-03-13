#include <iostream>
#include "MortgageBackedSecurities.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void RunQn1(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double time){

    double OAS = 0;
    bool isOAS = false;
    cout << "Please Enter the Kappa, rbar, and sigma" << endl;
    cin >> specialK >> rbar>> sigma;

    MortgageBackedSecurities::getNumerixPrepaymentModel( WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time, OAS, isOAS);

}

double RunQn2(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double time, double OAS, bool isOAS){

    double mbsPrice = MortgageBackedSecurities::getNumerixPrepaymentModel(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time, OAS, isOAS);
    cout << " Price is "<< mbsPrice << endl;
    return OAS;
}

void RunQn3(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double time, double OAS,
            bool isOAS){
    double p0 = 110000;
    double y = 0.0005;

    double p_plus = MortgageBackedSecurities::getNumerixPrepaymentModel( WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time, (OAS + y), isOAS);
    double p_minus = MortgageBackedSecurities::getNumerixPrepaymentModel( WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time, (OAS - y), isOAS);

    double duration = (p_minus - p_plus)/(2 * y * p0);
    double convexity = (p_plus + p_minus - 2 * p0)/(2 * p0 * y * y);

    cout <<"Duration: " << duration << endl;
    cout <<"Convexity: " << convexity << endl;
}


int main() {
    double WAC = 0.08;
    double pv0 = 100000;
    double r0 = 0.078;
    double specialK = 0.6;
    double rbar = 0.08;
    double sigma = 0.12;
    double time = 30; //num of years
    double OAS = 0;
    bool isOAS = true;

//    RunQn1(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time);

    OAS = -0.008492;
//    RunQn2(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time, OAS, isOAS);

    RunQn3(WAC, pv0, r0, specialK, rbar, sigma, time, OAS, isOAS);

    return 0;
}