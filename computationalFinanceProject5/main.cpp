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


void RunQn1(int nPath, int nSims, double sigma, double r, double x, double t, double s, double *tArray, double delta) {
    int k = 2; // k is set to 2 but can be 3 or 4

    int halfPath = nPath / 2;
    // matrix to show price path
    vector<vector<double> > priceProcess(nPath, vector<double>(nSims + 1, 0));
    // index matrix to understand when to exercise
    vector<vector<double> > index(nPath, vector<double>(nSims + 1, 0));
    // matrix to save the payoff
    vector<vector<double> > cashFlowMatrix(nPath, vector<double>(nSims + 1, 0));
    // matrix to save coefficent of the linear regression
    vector<vector<double> > a_coef(k, vector<double>(nSims, 0));

    //Generate Price process
    OptionPricing::generatePricePath(tArray, nSims, halfPath, s, r, sigma, priceProcess);

    // Determine the payoff and index and terminal stage
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

    // TODO: change the i to 0 for actual running
    for(int i=nSims - 1; i>=nSims - 1; i--){
        vector<vector< double > > A ( k, vector<double> ( k, 0 ) );
        vector<vector< double > > b ( k, vector<double> ( 1, 0 ) );
        vector<vector< double > > invMatrix ( k, vector<double> ( 1, 0 ) );
        // find out Matrix A for the current simulation Path
        OptionPricing::calculateMatrixA(priceProcess, A, k, nPath, i, 1);

        // find out Matrix b for the current simulation path
        OptionPricing::calcualateMatrixb(priceProcess, b, index, k, nPath, i, 1, nSims, r, delta, x);

        
        for(int i =0; i< k; i++){
            for(int j =0; j< k; j++){
                cout << A[i][j] <<",";
            }
            cout <<endl;
        }
    }


    cout << "#################################### price path ###################################" << endl;
    for (int i = 0; i < nPath; i++) {
        for (int j = 0; j <= nSims; j++) {
            cout << priceProcess[i][j] << ", ";
        }
        cout << endl;
    }

    cout << "#################################### cashflow matrix ###################################" << endl;
    for (int i = 0; i < nPath; i++) {
        for (int j = 0; j <= nSims; j++) {
            cout << cashFlowMatrix[i][j] << ", ";
        }
        cout << endl;
    }

    cout << "#################################### index ###################################" << endl;
    for (int i = 0; i < nPath; i++) {
        for (int j = 0; j <= nSims; j++) {
            cout << index[i][j] << ", ";
        }
        cout << endl;
    }


}


int main(){

    cout << "#################################### Qn 1 ###################################" << endl;
    int nPath = 10; // rows
    int nSims = 10; //cols
    double sigma = 0.2;
    double r = 0.06;
    double x = 40;
    double t = 0.5;
    double s = 40;
    double tArray[nSims];
    double delta = 1.0/nSims;
    // initialize to get delta t to generate weiner process
    for(int i=0; i<=nSims; i++){
        tArray[i] = delta * i * t;
    }

    RunQn1(nPath, nSims, sigma, r, x, t, s, tArray, delta);

    return 0;
}

