#include <iostream>
#include <cmath>
#include <algorithm>
#include "RandomGenerator.h"
#include "Mutils.h"

using namespace std;



void RunQn1(long double * z1, long double * z2, int size, double rho){

    long double * x;
    long double * y;
    double rho_sim;

    x = RandomGenerator::bivariateNormalX(z1, size);
    y = RandomGenerator::bivariateNormalY(z1, z2, rho, size);

    rho_sim = Mutils::Cov(x,y, size);
    cout << "#################################### Qn 1 ###################################" << endl;
    cout << "value of RHO using simulation is: " << rho_sim << endl;
    cout << "#############################################################################" << endl;
    cout << endl;
}

void RunQn2(long double * z1, long double * z2,  int size, double rho){
    long double * x;
    long double * y;

    x = RandomGenerator::bivariateNormalX(z1, size);
    y = RandomGenerator::bivariateNormalY(z1, z2, rho, size);

    long double arrE[size];

    for(int i = 0; i< size; i++){
        double val = double(x[i]*x[i]*x[i] + sin(y[i]) + x[i] * x[i] * y[i]);
        arrE[i] = max(0.0, val);
    }
    cout << "#################################### Qn 2 ###################################" << endl;
    cout << "Expected Value is : " << Mutils::Mean(arrE, size) << endl;
    cout << "#############################################################################" << endl;
    cout << endl;
}

void RunQn3(){

}

int main() {
    int size = 1000;
    double rho = -0.7;
    long double * uniArray = RandomGenerator::runif(size * 2);
    long double * z1;
    long double * z2;

    z1 = RandomGenerator::boxmuller(uniArray, size * 2);
    z2 = &z1[size];

    RunQn1(z1, z2, size, rho);
    rho = 0.6;
    RunQn2(z1, z2, size, rho);
    return 0;
}

