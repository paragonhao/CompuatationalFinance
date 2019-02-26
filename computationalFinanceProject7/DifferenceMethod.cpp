//
// Created by paragonhao on 25/2/19.
//

#include "DifferenceMethod.h"
#include "OptionPricing.h"
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include "Mutils.h"

using namespace Eigen;
using namespace std;

double DifferenceMethod::getPUEFD(double dt, double sigma, double delta_X, double r){
    return dt * ((sigma * sigma)/ (2 * delta_X * delta_X) + (r - 0.5 * sigma * sigma)/ (2 * delta_X));
}

double DifferenceMethod::getPMEFD(double dt, double sigma, double delta_X, double r){
    return 1 - (dt * sigma * sigma / (delta_X * delta_X)) - r * dt;
}

double DifferenceMethod::getPDEFD(double dt, double sigma, double delta_X, double r){
    return dt * ((sigma * sigma)/ (2 * delta_X * delta_X) - (r - 0.5 * sigma * sigma)/ (2 * delta_X));
}

void DifferenceMethod::EFDSolver(double currPrice, int deltaFactor){
    double s0 = currPrice;
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    double delta_x = sigma * sqrt(deltaFactor * dt);
    int time = int(t/dt) + 1; // cols
    int N = int(((log(30) -log(0.5))/delta_x)); // rows

    double pu = DifferenceMethod::getPUEFD(dt, sigma, delta_x, r);
    double pm = DifferenceMethod::getPMEFD(dt, sigma, delta_x, r);
    double pd = DifferenceMethod::getPDEFD(dt, sigma, delta_x, r);

    int totalPath = N * 2 + 1;

    //Initialize all needed matrix
    VectorXd logStockPath(totalPath);
    logStockPath = VectorXd::Zero(totalPath);

    VectorXd stockPrice(totalPath);
    stockPrice = VectorXd::Zero(totalPath);

    VectorXd payoffVectorF(totalPath);
    payoffVectorF = VectorXd::Zero(totalPath);

    VectorXd B(totalPath);
    B = VectorXd::Zero(totalPath);

    MatrixXd matA(totalPath, totalPath);
    matA = MatrixXd::Zero(totalPath, totalPath);

    // Generate Log Price path, price path and payoffVectorF,  total Path is the total number of rows.
    int counter = 0;
    for(int i = totalPath - 1; i>=0; i--){
        logStockPath(counter) = (i - N) * delta_x + log(s0);
        stockPrice(counter) =  exp(logStockPath(counter));
        payoffVectorF(counter) = Mutils::max(k - exp(logStockPath(counter)), 0) ;
        counter ++;
    }

    // Generate Matrix A

    //initialize the first and last row
    matA(0,0) = pu;
    matA(0,1) = pm;
    matA(0,2) = pd;
    matA(totalPath - 1,totalPath - 1) = pd;
    matA(totalPath - 1,totalPath - 2) = pm;
    matA(totalPath - 1,totalPath - 3) = pu;

    // start from the second row and end at second last row
    int startPos = 0;
    for(int i = 1; i < totalPath - 1; i++) {
        matA(i, startPos) = pd;
        matA(i, startPos + 1) = pm;
        matA(i, startPos + 2) = pu;
        startPos++;
    }


    B(totalPath-1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));

    // from the second and the second last
    for(int i = time; i > 0; i--){
        // New pay off
        payoffVectorF = matA * payoffVectorF + B ;
    }

    int size = int(payoffVectorF.size());

    int midpoint = (deltaFactor == 1)?(size/2 - 2):(size/2 - 1);

    for(int i =4; i<=16; i++){
        int ceilPos = midpoint - ceil(log(i/s0)/delta_x);
        int floorPos = midpoint - floor(log(i/s0)/delta_x);
        cout << i <<", " <<OptionPricing::putOptionPriceBS(r, sigma, t, i, k) << ", "<< (payoffVectorF(ceilPos) + payoffVectorF(floorPos)) / 2 << endl;
    }
}