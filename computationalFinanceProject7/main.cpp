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

void RunQn1(double currPrice){
    double s0 =currPrice;
    double sigma = 0.2;
    double k = 10;
    double r =0.04;
    double dt =0.002;
    double t = 0.5;

    double delta_x = sigma * sqrt(dt);
    int time = int(t/dt);
    int N = int(((log(17) -log(1))/delta_x));

    double pu = DifferenceMethod::getPUEFD(dt, sigma, delta_x, r);
    double pm = DifferenceMethod::getPMEFD(dt, sigma, delta_x, r);
    double pd = DifferenceMethod::getPDEFD(dt, sigma, delta_x, r);

    int totalPath = N * 2 + 1;
    VectorXd logStockPath(totalPath);
    VectorXd stockPrice(totalPath);
    VectorXd payoffVectorF(totalPath);

    // Generate Log Price path, price path and payoffVectorF,  total Path is the total number of rows.
    int counter = 0;
    for(int i = totalPath - 1; i>=0; i--){
        logStockPath(counter) = (i - N) * delta_x + log(s0);
        stockPrice(counter) =  exp(logStockPath(counter));
        payoffVectorF(counter) = Mutils::max(k - exp(logStockPath(counter)), 0) ;
        counter ++;
    }

    // Generate Matrix A
    MatrixXd matA(totalPath, totalPath);
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

    VectorXd B(totalPath);

    B(totalPath-1) = -(stockPrice(totalPath - 1) - stockPrice(totalPath - 2));

    // from the second and the second last
    for(int i = time; i > 0; i--){
        // New pay off
        payoffVectorF = matA * payoffVectorF + B ;
    }

    int size = int(payoffVectorF.size());
    cout << " totalPath " << totalPath << endl;
    cout <<payoffVectorF(size/2 - 1)<< endl;
//    for(int i =0; i< size; i++){
//        cout << i<<" "<<payoffVectorF(i)<<endl;
//    }


}


int main() {
    cout << "######################################## Qn1 ##################################################" <<endl;
    RunQn1(8);
    cout << "###############################################################################################" <<endl;
    return 0;
}