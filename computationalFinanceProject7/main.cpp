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

    int totalPath = int (currPrice * 2 / ds) + 1; // rows
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
            payoffVectorF(counter) = Mutils::max(stockPrice(counter) -k , 0);
            counter++;
        }

    }else if(type == "Put"){
        int counter = 0;
        for(int j= totalPath-1; j>=0;j--){
            stockPrice(counter) = j * ds;
            payoffVectorF(counter) = Mutils::max(k - stockPrice(counter) , 0);
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

    // initialize all as and bs, if price is 10, totalPath is 40
    // vector of as and bs are one less than the totalPath because we dont need it at maturity
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
    b3 = VectorXd::Zero(totalPath - 1);

    // as and bs will gointo matrix A
    int counter = 0;
    //assume currPrice = 10, totalPath = 41 (0-40),
    // we want to start from 39, in total 40 iterations (0-39)
    for(int j=totalPath - 2; j>=0; j--){

        double sigmaPow2JPow2 = sigma * sigma * j * j;

        a1(counter) = 0.5 * (sigmaPow2JPow2 - r * j) * (1 - alpha);
        a2(counter) = - (1/dt) - (sigmaPow2JPow2 + r) * (1 - alpha);
        a3(counter) = 0.5 * (sigmaPow2JPow2 + r * j) * (1 - alpha);

        b1(counter) = 0.5 * (sigmaPow2JPow2 - r * j) * alpha;
        b2(counter) = (1/dt) - (sigmaPow2JPow2 + r) * alpha;
        b3(counter) = 0.5 * (sigmaPow2JPow2 + r * j) * alpha;
        counter++;
    }


    // initialize MatrixA
    // totalPath= 41, 41 * 41 Matrix
    // last row last col: matA(totalPath - 1,totalPath - 1)
    MatrixXd matA(totalPath, totalPath);
    matA = MatrixXd::Zero(totalPath, totalPath);
    matA(0,0) = 1;
    matA(0,1) = -1;

    matA(totalPath - 1, totalPath - 1) = 1;
    matA(totalPath - 1, totalPath - 2) = -1;

    int startPos=0;
    for(int i = 1; i< totalPath - 1; i++){
        matA(i, startPos) = a1(i -1);
        matA(i, startPos + 1) = a2(i -1);
        matA(i, startPos + 2) = a3(i -1);
        startPos++;
    }

    // Inverse of A
    MatrixXd matAInverse(totalPath, totalPath);
    matAInverse = MatrixXd::Zero(totalPath, totalPath);
    matAInverse = matA.inverse();

    // initialize Matrix B
    MatrixXd matB(totalPath, totalPath);
    matB = MatrixXd::Zero(totalPath, totalPath);

    // Bs will going to the matrix
    int matBCounter = 0;
    for(int i=1; i < totalPath - 1; i++){
        matB(i, matBCounter) = -b1(i - 1);
        matB(i, matBCounter + 1) = -b2(i - 1);
        matB(i, matBCounter + 2) = -b3(i - 1);
        matBCounter++ ;
    }

    VectorXd vectorD(totalPath);
    vectorD = VectorXd::Zero(totalPath);

    // initialize the RHS of the equation
    vectorD = matB * payoffVectorF;
    // Initialize vectorD
    if(type == "Call"){
        vectorD(0) = stockPrice(0) - stockPrice(1);
        vectorD(totalPath - 1) = 0;
    }else if(type == "Put"){
        vectorD(0) = 0;
        vectorD(totalPath - 1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
    }
//    cout << vectorD << endl;

    for(int i = time; i > 0; i--){
        // getting pay off
        payoffVectorF = VectorXd::Zero(totalPath);
        payoffVectorF = matAInverse * vectorD;

        // American option, get pay off and compare with terminal pay off
        for(int j =0; j< totalPath; j++){
            payoffVectorF(j) = max(payoffVectorF(j),termialPayOff(j));
        }
        vectorD = matB * payoffVectorF;
        if(type == "Call"){
            vectorD(0) = stockPrice(0) - stockPrice(1);
            vectorD(totalPath - 1) = 0;
        }else if(type == "Put"){
            vectorD(0) = 0;
            vectorD(totalPath - 1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));
        }
    }

    int size = int(payoffVectorF.size());

    int midpoint = (size - 1)/2;

    cout << payoffVectorF(midpoint)<<endl;
}

int main() {

//    RunQn1EFD();
//    RunQn1IFD();
//    RunQn1CNFD();

    double ds = 0.25;
    string method = "CNFD";
//    string type = "Call";

    for(int i =4; i<=16; i ++ ){
        DifferenceMethod::GeneralisationOptionPriceSolver(i, ds, method);
    }


    return 0;
}