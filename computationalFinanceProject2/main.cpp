#include <iostream>
#include <cmath>
#include <algorithm>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"
#include <array>
#include <fstream>

void RunQn5b(int seed, double sigma, double s0, double r, string filename);

using namespace std;


double derivationVal(double x){
    return sqrt(1 - x * x);
}

long double t_x_function(long double x){
    return (1 - 0.74 * x * x )/(1- 0.74/3);
}

void RunQn1(long double * z1, long double * z2, int size, double rho){

    long double * x;
    long double * y;
    double rho_sim;

    x = RandomGenerator::bivariateNormalX(z1, size);
    y = RandomGenerator::bivariateNormalY(z1, z2, rho, size);

    rho_sim = Mutils::Corr(x, y, size);
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

void controlVar(long double * w_t, int size, double t ,int order){

    long double x[size];

    // x is used for EX, y is used for EY, w_t is the wiener process
    std::copy(w_t, w_t + size, x);
    long double y[size]; // y is cos(w)

    for(int i=0; i< size; i++){
        x[i] = cos(w_t[i]) * exp(t/2.0);
        y[i] = w_t[i] * w_t[i];
    }

    double gamma_t5 = Mutils::Cov(x, y, size) / Mutils::Cov(y, y, size);

    long double tArr[size];

    for(int i= 0; i< size; i++) {
        tArr[i] = x[i] - gamma_t5 * (y[i] - t);
    }
    cout << "Ea"<< order<<" mean before variance reduction: " << Mutils::Mean(x, size) << endl;
    cout << "Ea"<< order<<" sd before variance reduction: " << Mutils::StDev(x, size) << endl;

    cout << "Ea"<< order<<" mean after variance reduction: " << Mutils::Mean(tArr, size) << endl;
    cout << "Ea"<< order<<" sd after variance reduction: " << Mutils::StDev(tArr, size) << endl;
    cout << "#############################################################################" << endl;

}

void RunQn3(int seed){
    cout <<"#################################### Qn 3 (a) & (c) ###################################"<< endl;
    // first expectation
    double t = 5;
    double delta =0.001;

    int size = int(t/delta);
    long double * w_t = RandomGenerator::wienerProcess(t, size, seed);
    long double w_t5[size];

    // w_t5 is used for EX, y_t5 is used for EY
    std::copy(w_t, w_t + size, w_t5);
    long double y_t5[size];

    for(int i =0; i< size; i++){
        w_t5[i] = w_t[i] * w_t[i] + sin(w_t[i]);
        y_t5[i] = w_t[i] * w_t[i];
    }
    double t_5mean = Mutils::Mean(w_t5, size);
    double t_5sd = Mutils::StDev(w_t5,size);

    double gamma_t5 = Mutils::Cov(w_t5, y_t5, size) / Mutils::Cov(y_t5, y_t5, size);

    long double t5Arr[size];

    for(int i= 0; i< size; i++) {
        t5Arr[i] = w_t5[i] - gamma_t5 * (y_t5[i] - 5);
    }
    cout << "Ea1 mean before variance reduction: " << t_5mean << endl;
    cout << "Ea1 sd before variance reduction: " <<t_5sd << endl;

    cout << "Ea1 mean after variance reduction: " << Mutils::Mean(t5Arr, size) << endl;
    cout << "Ea1 sd before variance reduction: " <<Mutils::StDev(t5Arr, size) << endl;
    cout << "#############################################################################" << endl;


    // second expectation
    t = 0.5;
    size = int(t/delta);
    long double * w_t_1 = RandomGenerator::wienerProcess(t, size, seed+1);
    controlVar(w_t_1, size, t, 2);



    t = 3.2;
    size = int(t/delta);
    long double * w_t_2 = RandomGenerator::wienerProcess(t, size, seed+10);
    controlVar(w_t_2, size, t, 3);


    t = 6.5;
    delta =0.0001;
    size = int(t/delta);
    long double * w_t_3 = RandomGenerator::wienerProcess(t, size, seed+20);
    controlVar(w_t_3, size, t, 4);

    cout <<"#################################### Qn 3 (b) ###################################"<< endl;
    cout << "Expected Value at t=0.5, 3.2, 6.5 are close to 1" << endl;
    cout << "################################################################################" << endl;
}

void RunQn4(int seed){
    double r = 0.04;
    double sigma = 0.2;
    double s0 = 88;
    double T = 5;
    int size = 10000;
    double x = 100;

    cout <<"#################################### Qn 4(a) ###################################"<< endl;
    long double *w_t = RandomGenerator::wienerProcess(T,size,seed);

    long double *s_T = OptionPricing::callOptionPriceSimulation(w_t, size, r, sigma, T, s0, x);
    cout << "Call Option Price is: " << Mutils::Mean(s_T, size) << endl;

    cout <<"#################################### Qn 4(b) ###################################"<< endl;
    cout << "Call Option Price Computed by Black shole formula: " <<  OptionPricing::callOptionPriceBS(r, sigma, T, s0, x) << endl;


    // Using Reduction technique Antithetic Variable for monte carlo simulation
    long double y[size];

    std::copy(w_t, w_t + size, y);
    Mutils::Mutiply(y, size, -1); // we know that the variance is 5 in this case and expectaion value is 0

    long double *s_T_2 = OptionPricing::callOptionPriceSimulation(y, size, r, sigma, T, s0, x);

    long double tArr[size];

    for(int i=0; i<size; i++){
        tArr[i] = 0.5 * (s_T[i] + s_T_2[i]);
    }
    cout << "Call Option Price is: " << Mutils::Mean(tArr, size)<< endl;

}

void RunQn5a(int seed){
    cout <<"#################################### Qn 5(a) ###################################"<< endl;

    int size =10000, n=11;
    double r = 0.04;
    double sigma = 0.18;
    double s0 = 88;
    long double ESn[n];
    ESn[0] = s0;

    for(int i=0; i< n; i++){
        //w_i generation
        long double *arr = RandomGenerator::wienerProcess(i + 1, size, seed);
        double sum = 0;

        for(int j=0; j< size; j++){
           sum += s0 * exp(sigma * arr[j] + (r - sigma * sigma * 0.5) * (i + 1));
        }

        ESn[i+1] = sum / size ;
    }


    cout <<"################################# Write data to CSV file ...#####################"<< endl;
    Mutils::WriteToCSV(ESn, n, "../Data/Q5a1.csv");
    cout <<"################################################################################"<< endl;

}

void RunQn5b(int seed, double sigma, double s0, double r, string filename){
    cout <<"################################### Qn 5b #########################################"<< endl;
    int row = 6, col = 1001, pathNum = 1000;
    double delta = 10.0/1000;
    int size_5_b = row * col;
    long double s[row][col];
    long double * stdNor = RandomGenerator::boxmuller(RandomGenerator::runif(size_5_b, seed), size_5_b);
    int stdNorCounter = 0;

    for(int i = 0; i< row; i++){
        s[i][0] = s0;
        for(int j = 0; j < pathNum; j++){
            s[i][j+1] =  s[i][j] * exp( sigma * stdNor[stdNorCounter] * sqrt(delta)  +  r * delta);
            stdNorCounter ++;
        }
    }

    // Writing to CSV
    ofstream file;
    file.open(filename);
    for(int i=0; i<row; i++){
       for(int j=0; j<col; j++){
            file << s[i][j] << ",";
       }
       file << "\n";
    }
    cout <<"##############################################################################"<< endl;
}

void RunQn5d(int seed){

    int size =10000, n=11;
    double r = 0.04;
    double sigma = 0.35;
    double s0 = 88;
    long double ESn[n];
    ESn[0] = s0;

    for(int i=0; i< n; i++){
        //w_i generation
        long double *arr = RandomGenerator::wienerProcess(i + 1, size, seed);

        double sum = 0;
        for(int j=0; j< size; j++){
            sum += s0 * exp(sigma * arr[j] + (r - sigma * sigma * 0.5) * (i + 1));
        }

        ESn[i+1] = sum / size ;
    }
    cout <<"################################# Write data to CSV file ...#####################"<< endl;
    Mutils::WriteToCSV(ESn, n, "../Data/Q5dESN.csv");
    cout <<"################################################################################"<< endl;
}

void RunQn6a(){
    // use euler method
    double size = 1000;
    double h = 0.001;
    //initial value
    double x0 = 0.0;
    double y0= 1.0;
    double y;

    for(int i =0;i <size; i++){
        y = y0 + h * derivationVal(x0);
        y0 = y;
        x0 += 0.001;
    }
    cout <<"#################################### Qn 6(a) ###################################"<< endl;
    cout << "Euler method to get the integral, value is: "<< (y -1) * 4 << endl;
    cout <<"################################################################################"<< endl;
}

void RunQn6b(int seed){
    int size = 10000;
    long double * uniArray = RandomGenerator::runif(size, 1234);
    long double gArr[size];

    for(int i=0; i< size; i++){
        gArr[i] = sqrt(1 - uniArray[i] * uniArray[i]);
    }
    cout <<"#################################### Qn 6(b) ###################################"<< endl;
    cout << "Monte Carlo Approx is: "<< Mutils::Mean(gArr,size) * 4 << endl;
    cout <<"################################################################################"<< endl;
}

void RunQn6c(int seed){
    // generate t(x) distribution
    long double *yArr;
    int counter = 0;
    long double y_rv =0;
    long double sum =0.0;
    int size =10000;
    long double *y = RandomGenerator::runif(size, seed);
    long double *uni = RandomGenerator::runif(size, seed + 1000);
    long double result = 0.0;

    for(int i=0; i< size; i++){
        y_rv = t_x_function(y[i]);
        if(uni[i] <= y_rv/1.5){
            result +=derivationVal(uni[i])/t_x_function(uni[i]);
            counter++;
        }
    }
    cout <<"#################################### Qn 6(b) ###################################"<< endl;
    cout << "Monte Carlo Using Importance Sampling Method: "<< result/counter *4<< endl;
    cout <<"################################################################################"<< endl;
}

int main() {
    int size = 1000;
    double rho = -0.7;
    const int seed = 1234567890;
    long double * uniArray = RandomGenerator::runif(size * 2, seed);
    double r = 0.04;
    double sigma = 0.18;
    double s0 = 88;

    long double * z1;
    long double * z2;

    z1 = RandomGenerator::boxmuller(uniArray, size * 2);
    z2 = &z1[size];

    RunQn1(z1, z2, size, rho);
    rho = 0.6;
    RunQn2(z1, z2, size, rho);
    RunQn3(seed);
    RunQn4(seed);
    RunQn5a(seed);
    RunQn5b(seed, sigma, s0, r,"../Data/Q5b1.csv");
    sigma = 0.35;
    RunQn5b(seed, sigma, s0, r,"../Data/Q5dhigsigma.csv");
    RunQn5d(seed);
    RunQn6a();
    RunQn6b(seed);
    RunQn6c(seed);
    return 0;
}

