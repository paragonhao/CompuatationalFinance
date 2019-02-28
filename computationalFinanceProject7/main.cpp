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

void RunQn2(){

    cout << "######################################## Qn2 ##################################################" <<endl;
    cout << "ds = 0.25"<< endl;
    double ds = 0.25;
    for(int i =4; i<=16; i ++ ){
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "EFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "IFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "CNFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "EFD", "Put");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "IFD", "Put");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "CNFD", "Put");
    }


    cout << "ds = 1"<< endl;
    ds =1;
    for(int i =4; i<=16; i ++ ){
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "EFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "IFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "CNFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "EFD", "Put");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "IFD", "Put");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "CNFD", "Put");
    }

    cout << "ds = 1.25"<< endl;
    ds =1.25;
    for(int i =4; i<=16; i ++ ){
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "EFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "IFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "CNFD", "Call");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "EFD", "Put");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "IFD", "Put");
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "CNFD", "Put");
    }

}

int main() {

//    RunQn1EFD();
//    RunQn1IFD();
//    RunQn1CNFD();
//    RunQn2();
    double ds =0.25;
    for(int i =4; i<=16; i ++ ){
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, "CNFD", "Call");
    }

    return 0;
}