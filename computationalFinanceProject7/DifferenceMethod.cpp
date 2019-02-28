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

// EFD
double DifferenceMethod::getPUEFD(double dt, double sigma, double delta_X, double r){
    return dt * ((sigma * sigma)/ (2 * delta_X * delta_X) + (r - 0.5 * sigma * sigma)/ (2 * delta_X));
}

double DifferenceMethod::getPMEFD(double dt, double sigma, double delta_X, double r){
    return 1 - (dt * sigma * sigma / (delta_X * delta_X)) - r * dt;
}

double DifferenceMethod::getPDEFD(double dt, double sigma, double delta_X, double r){
    return dt * ((sigma * sigma)/ (2 * delta_X * delta_X) - (r - 0.5 * sigma * sigma)/ (2 * delta_X));
}

// IFD
double DifferenceMethod::getPUIFD(double dt, double sigma, double delta_X, double r){
    return -0.5 * dt * ((sigma*sigma)/(delta_X*delta_X) + (r - 0.5 *sigma*sigma)/delta_X );
}

double DifferenceMethod::getPMIFD(double dt, double sigma, double delta_X, double r){
    return 1 + dt * ((sigma * sigma)/(delta_X * delta_X)) + r * dt;
}

double DifferenceMethod::getPDIFD(double dt, double sigma, double delta_X, double r){
    return -0.5 * dt * ((sigma*sigma)/(delta_X*delta_X) - (r - 0.5 *sigma*sigma)/delta_X);
}

//CNFD
double DifferenceMethod::getPUCNFD(double dt, double sigma, double delta_X, double r){
    return -0.25* dt * ((sigma*sigma)/(delta_X*delta_X) + (r - 0.5* sigma * sigma)/delta_X);
}

double DifferenceMethod::getPMCNFD(double dt, double sigma, double delta_X, double r){
    return 1 + dt * sigma * sigma * 0.5 / (delta_X * delta_X) + r * dt * 0.5;
}

double DifferenceMethod::getPDCNFD(double dt, double sigma, double delta_X, double r){
    return -0.25* dt * ((sigma*sigma)/(delta_X*delta_X) - (r - 0.5* sigma * sigma)/delta_X);
}



void DifferenceMethod::EFDEuroPutSolver(double currPrice, int deltaFactor){
    double s0 = currPrice;
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    double delta_x = sigma * sqrt(deltaFactor * dt);
    int time = int(t/dt) + 1; // cols
//    int N = int(((log(currPrice * 3 * sigma + currPrice) -log(currPrice - currPrice * 3 * sigma))/delta_x)); // rows
    int N = 50;
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
        matA(i, startPos) = pu;
        matA(i, startPos + 1) = pm;
        matA(i, startPos + 2) = pd;
        startPos++;
    }


    B(totalPath-1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));

    // from the second and the second last
    for(int i = time; i > 0; i--){
        // New pay off
        payoffVectorF = matA * payoffVectorF + B ;
    }

    int size = int(payoffVectorF.size());

    int midpoint = (size - 1)/2;

    //price for mid point
    cout << currPrice<<", "<<payoffVectorF(midpoint)<<", "<<OptionPricing::putOptionPriceBS(r, sigma, t, currPrice, k) <<endl;

}

void DifferenceMethod::IFDEuroPutSolver(double currPrice, int deltaFactor){
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

    // Before running the loop, initialize the inverese A
    MatrixXd matAInverse(totalPath, totalPath);
    matAInverse = MatrixXd::Zero(totalPath, totalPath);
    matAInverse = matA.inverse();

    for(int i = time; i > 0; i--){
        payoffVectorF = matAInverse * B;
        B = payoffVectorF;
        B(totalPath-1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
        B(0) = 0;
    }

    int size = int(payoffVectorF.size());

    int midpoint = (size - 1)/2;

    cout << currPrice<<", "<<payoffVectorF(midpoint)<<", "<<OptionPricing::putOptionPriceBS(r, sigma, t, currPrice, k) <<endl;
}

void DifferenceMethod::CNFDEuroPutSolver(double currPrice, int deltaFactor){
    double s0 = currPrice;
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    double delta_x = sigma * sqrt(deltaFactor * dt);
    int time = int(t/dt) + 1; // cols

    int N = 80;

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


    // Before running the loop, initialize the inverese A
    MatrixXd matAInverse(totalPath, totalPath);
    matAInverse = MatrixXd::Zero(totalPath, totalPath);
    matAInverse = matA.inverse();

//    cout << vectorZ<<endl;
    for(int i = time; i > 0; i--){
        // next pay off
        payoffVectorF = matAInverse * vectorZ;
        // create new vector Z with the current pay off;
        vectorZ = matZMutilply * payoffVectorF;
        vectorZ(totalPath - 1)  = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
    }

    int size = int(payoffVectorF.size());

    int midpoint = (size - 1)/2;

    cout << currPrice<<", "<<payoffVectorF(midpoint)<<", "<<OptionPricing::putOptionPriceBS(r, sigma, t, currPrice, k) <<endl;
}

// Generalisation Solver
void DifferenceMethod::GeneralisationOptionPriceSolver(double currPrice, double ds, string method, string type){
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    int N = int(currPrice / ds);
    int totalPath = N * 2 + 1; // rows
    int time = int(t/dt) + 1; // cols

    VectorXd stockPrice(totalPath);
    stockPrice = VectorXd::Zero(totalPath);

    VectorXd payoffVectorF(totalPath);
    payoffVectorF = VectorXd::Zero(totalPath);

    VectorXd termialPayOff(totalPath);
    termialPayOff = VectorXd::Zero(totalPath);

    //######################## decide call put pay off ##########################

    if(type == "Call"){
        int counter = 0;
        for(int j= totalPath-1; j>=0;j--){
            stockPrice(counter) = j * ds;
            payoffVectorF(counter) = Mutils::max(stockPrice(counter) -k , 0.0);
            counter++;
        }

    }else if(type == "Put"){
        int counter = 0;
        for(int j= totalPath-1; j>=0;j--){
            stockPrice(counter) = j * ds;
            payoffVectorF(counter) = Mutils::max(k - stockPrice(counter) , 0.0);
            counter++;
        }
    }

    // make a copy of the terminal pay off
    termialPayOff = payoffVectorF;

    //############################################################################
    double alpha = 0;

    // choose the alpha according to the method
    if(method == "EFD"){
        alpha = 1;
    }else if(method == "IFD"){
        alpha = 0;
    }else if(method == "CNFD"){
        alpha = 0.5;
    }

    // initialize MatrixA
    // totalPath= 41, 41 * 41 Matrix
    // last row last col: matA(totalPath - 1,totalPath - 1)
    MatrixXd matA(totalPath, totalPath);
    matA = MatrixXd::Zero(totalPath, totalPath);
    matA(0,0) = 1;
    matA(0,1) = -1;

    matA(totalPath - 1, totalPath - 2) = 1;
    matA(totalPath - 1, totalPath - 1) = -1;

    MatrixXd matB(totalPath, totalPath);
    matB = MatrixXd::Zero(totalPath, totalPath);
    matB(0,0) = 1;
    matB(0,1) = -1;
    matB(totalPath - 1, totalPath - 2) = 1;
    matB(totalPath - 1, totalPath - 1) = -1;

    VectorXd a1(totalPath - 1);
    a1 = VectorXd::Zero(totalPath - 1);

    VectorXd a2(totalPath - 1);
    a2 = VectorXd::Zero(totalPath - 1);

    VectorXd a3(totalPath - 1);
    a3 = VectorXd::Zero(totalPath - 1);

    VectorXd b1(totalPath - 1);
    b1 = VectorXd::Zero(totalPath - 1);

    VectorXd b2(totalPath - 1);
    b2 = VectorXd::Zero(totalPath - 1);

    VectorXd b3(totalPath - 1);
    b3 = VectorXd::Zero(totalPath);

    int Acounter = 0;
    //assume currPrice = 10, totalPath = 41 (0-40),
    // we want to start from 39, in total 40 iterations (0-39)
    for(int i = 1; i< totalPath - 1; i++){

        double j = totalPath - 1 - i;
        double sigmaPow2JPow2 = sigma * sigma * j * j;

        a1(Acounter) = 0.5 * (sigmaPow2JPow2 - r * j) * (1 - alpha);
        a2(Acounter) = - (1/dt) - (sigmaPow2JPow2 + r) * (1 - alpha);
        a3(Acounter) = 0.5 * (sigmaPow2JPow2 + r * j) * (1 - alpha);

        b1(Acounter) = 0.5 * (sigmaPow2JPow2 - r * j) * alpha;
        b2(Acounter) = (1/dt) - (sigmaPow2JPow2 + r) * alpha;
        b3(Acounter) = 0.5 * (sigmaPow2JPow2 + r * j) * alpha;
        Acounter++;
    }

    int startPos=0;
    for(int i = 1; i< totalPath - 1; i++){
        matA(i, startPos) = a3(i -1);
        matA(i, startPos + 1) = a2(i -1);
        matA(i, startPos + 2) = a1(i -1);

        matB(i, startPos) = -b3(i - 1);
        matB(i, startPos + 1) = -b2(i - 1);
        matB(i, startPos + 2) = -b1(i - 1);
        startPos++;
    }

    // Inverse of A
    MatrixXd matAInverse(totalPath, totalPath);
    matAInverse = MatrixXd::Zero(totalPath, totalPath);
    matAInverse = matA.inverse();

    // initialize Matrix B
    VectorXd vectorD(totalPath);
    vectorD = VectorXd::Zero(totalPath);

    // initialize the RHS of the equation
    vectorD = matB * payoffVectorF;

    for(int i = time; i > 0; i--){

        if(type == "Call"){
            vectorD(0) = stockPrice(0) - stockPrice(1);
            vectorD(totalPath - 1) = 0;
            payoffVectorF = matAInverse * vectorD;
            // American option, get pay off and compare with terminal pay off
            for(int j =0; j< totalPath; j++){
                payoffVectorF(j) = max(payoffVectorF(j),termialPayOff(j));
            }
        }else if(type == "Put"){
            vectorD(0) = 0;
            vectorD(totalPath - 1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
            payoffVectorF = matAInverse * vectorD;
            // American option, get pay off and compare with terminal pay off
            for(int j =0; j< totalPath; j++){
                payoffVectorF(j) = max(payoffVectorF(j),termialPayOff(j));
            }
        }

        vectorD = matB * payoffVectorF;
    }

    int size = int(payoffVectorF.size());

    int midpoint = (size - 1)/2;

    cout << currPrice << ", "<< payoffVectorF(midpoint)<<endl;
}