#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"
#include "JDDefaultOption.h"
#include "DifferenceMethod.h"
#include <random>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void RunQn1EFD(){
    cout << "######################################## Qn1 ##################################################" <<endl;
    cout << "Explicit Finite-Difference method: sigma * sqrt(dx)" <<endl;

    int deltaFactor = 1;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    for(int i =4; i<=16; i++){
        DifferenceMethod::EFDEuroPutSolver(i, deltaFactor);
    }

    cout << "sigma * sqrt(3 * dx)" <<endl;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    deltaFactor = 3;
    for(int i =4; i<=16; i++){
        DifferenceMethod::EFDEuroPutSolver(i, deltaFactor);
    }

    cout << "sigma * sqrt(4 * dx)" <<endl;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    deltaFactor = 4;
    for(int i =4; i<=16; i++){
        DifferenceMethod::EFDEuroPutSolver(i, deltaFactor);
    }
    cout << "###############################################################################################" <<endl;
}

void RunQn1IFD(){
    cout << "######################################## Qn1 ##################################################" <<endl;
    cout << "Implicit Finite-Difference method: sigma * sqrt(dx)" <<endl;

    int deltaFactor = 1;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    for(int i =4; i<=16; i++){
        DifferenceMethod::IFDEuroPutSolver(i, deltaFactor);
    }

    cout << "sigma * sqrt(3 * dx)" <<endl;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    deltaFactor = 3;
    for(int i =4; i<=16; i++){
        DifferenceMethod::IFDEuroPutSolver(i, deltaFactor);
    }

    cout << "sigma * sqrt(4 * dx)" <<endl;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    deltaFactor = 4;
    for(int i =4; i<=16; i++){
        DifferenceMethod::IFDEuroPutSolver(i, deltaFactor);
    }

    cout << "###############################################################################################" <<endl;
}

void RunQn1CNFD(){
    cout << "######################################## Qn1 ##################################################" <<endl;
    cout << "Crank-Nicolson Finite Difference method: sigma * sqrt(dx)" <<endl;

    int deltaFactor = 1;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    for(int i =4; i<=16; i++){
        DifferenceMethod::CNFDEuroPutSolver(i, deltaFactor);
    }

    cout << "sigma * sqrt(3 * dx)" <<endl;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    deltaFactor = 3;
    for(int i =4; i<=16; i++){
        DifferenceMethod::CNFDEuroPutSolver(i, deltaFactor);
    }

    cout << "sigma * sqrt(4 * dx)" <<endl;
    cout << "Price"<<", "<<"Pay off"<<", BS"<<endl;
    deltaFactor = 4;
    for(int i =4; i<=16; i++){
        DifferenceMethod::CNFDEuroPutSolver(i, deltaFactor);
    }

    cout << "###############################################################################################" <<endl;
}

void RunQn2EFD(double currPrice, double ds, string method, string type){
    double s0 = currPrice;
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    int totalPath = int (currPrice * 2 / ds);
    cout << totalPath << endl;

    VectorXd stockPrice(totalPath);
    stockPrice = VectorXd::Zero(totalPath);



}

int main() {

//    RunQn1EFD();
//    RunQn1IFD();
//    RunQn1CNFD();


    double currPrice = 10;
    double ds = 0.5;
    string method = "EFD";
    string type = "Call";
    RunQn2EFD(currPrice, ds, method,type);


    return 0;
}