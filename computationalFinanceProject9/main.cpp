#include <iostream>
#include "MortgageBackedSecurities.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void RunQn1(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double duration){

    double OAS = 0;
    bool isOAS = false;
    cout << "Please Enter the Kappa, rbar, and sigma" << endl;
    cin >> specialK >> rbar>> sigma;

    MortgageBackedSecurities::getNumerixPrepaymentModel(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  duration, OAS, isOAS);

}

double RunQn2(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double duration, double OAS, bool isOAS){

    double mbsPrice = MortgageBackedSecurities::getNumerixPrepaymentModel(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  duration, OAS, isOAS);
    cout << " Price is "<< mbsPrice << endl;
    return OAS;
}

void RunQn3(double WAC, double pv0, double r0, double specialK, double rbar, double sigma, double duration, double OAS,
            bool isOAS){
    double p0 = 110000;
    double y = 0.0005;

    double p_plus = MortgageBackedSecurities::getNumerixPrepaymentModel( WAC,  pv0,  r0,  specialK,  rbar,  sigma,  duration, (OAS + y), isOAS);
    double p_minus = MortgageBackedSecurities::getNumerixPrepaymentModel( WAC,  pv0,  r0,  specialK,  rbar,  sigma,  duration, (OAS - y), isOAS);

    double mbsduration = (p_minus - p_plus)/(2 * y * p0);
    double convexity = (p_plus + p_minus - 2 * p0)/(2 * p0 * y * y);

    cout <<"Duration: " << mbsduration << endl;
    cout <<"Convexity: " << convexity << endl;
}


int main() {
    double WAC = 0.08;
    double pv0 = 100000;
    double r0 = 0.078;
    double specialK = 0.6;
    double rbar = 0.08;
    double sigma = 0.1;
    double duration = 30.0; //num of years
    double OAS = 0;
    bool isOAS = true;


    RunQn1(WAC,  pv0,  r0,  specialK,  rbar,  sigma ,  duration);

    OAS = -0.012811;
    RunQn2(WAC,  pv0,  r0,  specialK,  rbar,  sigma,  duration, OAS, isOAS);

    RunQn3(WAC, pv0, r0, specialK, rbar, sigma, duration, OAS, isOAS);

//For loop to get the values for the graph

//    for(int i = 0; i< 11;i++){
//        RunQn1(WAC,  pv0,  r0,  specialK,  rbar,  sigma + i * 0.01,  duration);
//    }
    return 0;
}