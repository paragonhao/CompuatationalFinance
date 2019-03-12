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

void RunQn2(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double time, double OAS, bool isOAS){

    OAS = -0.0088; // starting position
    double mbsPrice = MortgageBackedSecurities::getNumerixPrepaymentModel(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time, OAS, isOAS);
    cout << " Price is "<< mbsPrice << endl;

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
    RunQn2(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  time, OAS, isOAS);


    return 0;
}