#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include <chrono>
#include "RandomGenerator.h"
#include "Mutils.h"
#include <Eigen/Dense>
#include "OptionPricing.h"

using namespace std;
using namespace std::chrono;
using namespace Eigen;


void RunQn1(int nPath, int nSims, double sigma, double r, double x, double s, double *tArray, double delta) {
    int k = 2; // k is set to 2 but can be 3 or 4
    int method =1; // method is set to 1 but can be 1 or 2 or 3;

    int halfPath = nPath / 2;

    // matrix to show price path
    vector<vector<double> > priceProcess(nPath, vector<double>(nSims + 1, 0));
    // index matrix to understand when to exercise
    vector<vector<double> > index(nPath, vector<double>(nSims + 1, 0));
    // matrix to save the payoff
    vector<vector<double> > cashFlowMatrix(nPath, vector<double>(nSims + 1, 0));

    //Generate Price process
    OptionPricing::generatePricePath(tArray, nSims, halfPath, s, r, sigma, priceProcess);

//#################################### code to verify price process ###################################"
//    double total =0;
//    for (int i = 0; i < nPath; i++) {
//
//        total += Mutils::max(x - priceProcess[i][nSims],0);
//
//    }
//    cout << total/nPath<<endl;
//########################################################################################################################
    //Determine the payoff and index and terminal stage

    for(int j = nSims; j >=0; j--) {
        for (int i = 0; i < nPath; i++) {
            double EV = x - priceProcess[i][j];

            // At the last step, ECV = 0 and payoff is just EV;
            if (j == nSims) {
                if (EV > 0) {
                    cashFlowMatrix[i][j] = EV;
                    index[i][j] = 1;
                }
            }
            cashFlowMatrix[i][j] = (EV > 0) ? EV : 0;
        }
    }
//
//    // TOOD: change i>6 later to i>0
    for(int i=nSims - 1; i>0; i--){

        MatrixXd matA(k, k);
        VectorXd matb(k);
        VectorXd result(k);
        matA = MatrixXd::Zero(k, k);
        matb = VectorXd::Zero(k);

        // find out Matrix A for the current simulation Path
        OptionPricing::calculateMatrixA(priceProcess, matA, k, nPath, i, 1);

        // find out Matrix b for the current simulation path
        OptionPricing::calcualateMatrixb(cashFlowMatrix,priceProcess, matb, index, k, nPath, i, 1, nSims, r, delta, x);
        result = matA.inverse() * matb;


        //so far so good
        OptionPricing::calculateContinuationValue(result, priceProcess, index, cashFlowMatrix, nSims, i, nPath, k, 1);
    }
//
    double payoff = OptionPricing::calculateFinalPayOff(cashFlowMatrix, index, nPath, nSims, r, delta);
    cout <<"Continuation Value is: " <<payoff << endl;

//    cout << "#################################### Price Process matrix ###################################" << endl;
//    for (int i = 0; i < nPath; i++) {
//        for (int j = 0; j <= nSims; j++) {
//            cout << priceProcess[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//    cout << "#################################### cashflow matrix ###################################" << endl;
//    for (int i = 0; i < nPath; i++) {
//        for (int j = 0; j <= nSims; j++) {
//            cout << cashFlowMatrix[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//    cout << "#################################### index ###################################" << endl;
//    for (int i = 0; i < nPath; i++) {
//        for (int j = 0; j <= nSims; j++) {
//            cout << index[i][j] << " ";
//        }
//        cout << endl;
//    }


}


int main(){

    cout << "#################################### Qn 1 ###################################" << endl;
    int nPath =1000; // rows
    int nSims = 1000; //cols
    double sigma = 0.2;
    double r = 0.06;
    double x = 40;
    double t = 0.5;
    double s = 40;
    auto * tArray = new double [nSims+1];
    double delta = t/nSims;
    // initialize to get delta t to generate weiner process
    for(int i=0; i<=nSims; i++){
        tArray[i] = delta * i;
    }

    //American Put
    RunQn1(nPath, nSims, sigma, r, x, s, tArray, delta);

    return 0;
}

