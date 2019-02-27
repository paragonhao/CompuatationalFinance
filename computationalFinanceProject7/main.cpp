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
    cout << "Price, Black sholes, Put Pay off" <<endl;
    int deltaFactor = 1;
    DifferenceMethod::EFDSolver(10, deltaFactor);
    cout << "sigma * sqrt(3 * dx)" <<endl;
    deltaFactor = 3;
    DifferenceMethod::EFDSolver(10, deltaFactor);
    cout << "sigma * sqrt(4 * dx)" <<endl;
    deltaFactor = 4;
    DifferenceMethod::EFDSolver(10, deltaFactor);
    cout << "###############################################################################################" <<endl;
}

void RunQn1IFD(double currPrice, int deltaFactor){
    cout << "######################################## Qn1 ##################################################" <<endl;
    cout << "Implicit Finite-Difference method: sigma * sqrt(dx)" <<endl;
    double s0 = currPrice;
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    double delta_x = sigma * sqrt(deltaFactor * dt);
    int time = int(t/dt) + 1;
//    int N = int(((log(currPrice * 3 * sigma + currPrice) -log(currPrice - currPrice * 3 * sigma))/delta_x)); // rows
    int N = 50;
    double pu = DifferenceMethod::getPUIFD(dt, sigma, delta_x, r);
    double pm = DifferenceMethod::getPMIFD(dt, sigma, delta_x, r);
    double pd = DifferenceMethod::getPDIFD(dt, sigma, delta_x, r);

    int totalPath = N * 2 + 1; // total number of rows

    VectorXd logStockPath(totalPath);
    logStockPath = VectorXd::Zero(totalPath);

    VectorXd stockPrice(totalPath);
    stockPrice = VectorXd::Zero(totalPath);

    VectorXd payoffVectorF(totalPath);
    payoffVectorF = VectorXd::Zero(totalPath);

    VectorXd B(totalPath);
    B = VectorXd::Zero(totalPath);

    //initialize the first and last row
    MatrixXd matA(totalPath, totalPath);
    matA = MatrixXd::Zero(totalPath, totalPath);


    // Generate price process
    int counter = 0;
    for(int i = totalPath - 1; i>=0; i--){
        logStockPath(counter) = (i - N) * delta_x + log(s0);
        stockPrice(counter) =  exp(logStockPath(counter));
        payoffVectorF(counter) = Mutils::max(k - exp(logStockPath(counter)), 0) ;
        counter ++;
    }

    matA(0,0) = 1;
    matA(0,1) = -1;
    matA(totalPath - 1,totalPath - 1) = -1;
    matA(totalPath - 1,totalPath - 2) = 1;


    int startPos = 0;
    for(int i = 1; i < totalPath - 1; i++) {
        matA(i, startPos) = pu;
        matA(i, startPos + 1) = pm;
        matA(i, startPos + 2) = pd;
        startPos++;
    }

    B = payoffVectorF;
    B(totalPath-1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
    B(0) = 0;

    for(int i = time; i > 0; i--){
        cout << i << endl;
        payoffVectorF = matA.inverse() * B;
        B = payoffVectorF;
        B(totalPath-1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
        B(0) = 0;
    }

    int size = int(payoffVectorF.size());

    int midpoint = (size - 1)/2;

    cout << payoffVectorF(midpoint)<<endl;

    cout << "###############################################################################################" <<endl;
}


void RunQn1CNFD(double currPrice, int deltaFactor){
    cout << "######################################## Qn1 ##################################################" <<endl;
    cout << "Crank-Nicolson Finite-Difference method: sigma * sqrt(dx)" <<endl;
    double s0 = currPrice;
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    double delta_x = sigma * sqrt(deltaFactor * dt);
    int time = int(t/dt) + 1; // cols

    int N = 50;

    double pu = DifferenceMethod::getPUCNFD(dt, sigma, delta_x, r);
    double pm = DifferenceMethod::getPMCNFD(dt, sigma, delta_x, r);
    double pd = DifferenceMethod::getPDCNFD(dt, sigma, delta_x, r);

    int totalPath = N * 2 + 1;

    VectorXd logStockPath(totalPath);
    logStockPath = VectorXd::Zero(totalPath);

    VectorXd stockPrice(totalPath);
    stockPrice = VectorXd::Zero(totalPath);

    VectorXd payoffVectorF(totalPath);
    payoffVectorF = VectorXd::Zero(totalPath);

    int counter = 0;
    for(int i = totalPath - 1; i>=0; i--){
        logStockPath(counter) = (i - N) * delta_x + log(s0);
        stockPrice(counter) =  exp(logStockPath(counter));
        payoffVectorF(counter) = Mutils::max(k - exp(logStockPath(counter)), 0) ;
        counter ++;
    }

    // Initialize Matrix A
    MatrixXd matA(totalPath, totalPath);
    matA = MatrixXd::Zero(totalPath, totalPath);
    matA(0,0) = 1;
    matA(0,1) = -1;
    matA(totalPath - 1,totalPath - 1) = -1;
    matA(totalPath - 1,totalPath - 2) = 1;

    int startPos = 0;
    for(int i = 1; i < totalPath - 1; i++) {
        matA(i, startPos) = pu;
        matA(i, startPos + 1) = pm;
        matA(i, startPos + 2) = pd;
        startPos++;
    }

    //Initialize the matrix to mutiply pay off to get Z
    MatrixXd matZMutilply;
    matZMutilply = MatrixXd::Zero(totalPath, totalPath);

    startPos = 0;
    for(int i = 1; i < totalPath - 1; i++) {
        matZMutilply(i, startPos) = -pu;
        matZMutilply(i, startPos + 1) = -(pm -2);
        matZMutilply(i, startPos + 2) = -pd;
        startPos++;
    }

//    // Initialize the the first vector Z using the pay off at maturity and Z matrix,
    VectorXd vectorZ(totalPath);
    vectorZ = VectorXd::Zero(totalPath);

    vectorZ = matZMutilply * payoffVectorF;
    vectorZ(totalPath - 1)  = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));

//    cout << vectorZ<<endl;
    for(int i = time; i > 0; i--){
        // next pay off
        payoffVectorF = matA.inverse() * vectorZ;
        // create new vector Z with the current pay off;
        vectorZ = matZMutilply * payoffVectorF;
        vectorZ(totalPath - 1)  = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
    }

    int size = int(payoffVectorF.size());

    int midpoint = (size - 1)/2;

    cout << payoffVectorF(midpoint)<<endl;

}

int main() {

//    RunQn1EFD();
//    RunQn1IFD(8, 4);
    RunQn1CNFD(8, 3);
    return 0;
}