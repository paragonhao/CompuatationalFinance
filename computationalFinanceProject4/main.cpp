#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"

using namespace std;


void RunQn1(){
    double r = 0.05;
    double sigma = 0.24;
    double s0 = 32;
    double k = 30;
    double t = 0.5;

    vector<int> nSize {10, 20, 40, 80, 100, 200, 500};
    vector<string> method {"a", "b", "c", "d"};
    cout << "Option price using Black Sholes: " << OptionPricing::callOptionPriceBS(r,sigma,t, s0, k) << endl;

    for(int j=0; j<4; j++){
        for(int i=0; i< 7; i++){
            cout << "Option Price for method ("<< method[j]<< "), n = "<< nSize[i] <<" is: " << OptionPricing::callOptionEuropeanBinomial(method[j], s0, k, r, sigma, t, nSize[i]) << endl;
        }
        cout << "####################################################################################" << endl;
    }
}

void RunQn2(){
    double rf = 0.02;
    double price = 1141.42;
    double k = (int(price)/10) * 10.0 * 1.1;
    int steps = 200;

    ifstream  data("../Data/GOOG_monthly_ret.csv");


    string line;
    int counter = 0;
    long double *ret = new long double[60];

    while(getline(data,line))
    {
        long double temp = stold(line);
        ret[counter] = temp;
        counter ++;
    }

    long double sigma = Mutils::StDev(ret, 60) * sqrt(12);
    cout << "Price is 58.8 on Yahoo Finance"<< endl;
    cout << "Volatility based on calculation is: "<< sigma << endl;
    cout << "Estimate the option price on Jan 2020: " <<OptionPricing::callOptionEuropeanBinomial("b", price, k, rf, sigma, 0.968, 200)<< endl;
    cout << "Estimate the option price with Volatility equal to 20.4%: " <<OptionPricing::callOptionEuropeanBinomial("b", price, k, rf, 0.204, 1, 200)<< endl;

}

void RunQn3(){
    double s0= 49;
    double k = 50;
    double r = 0.03;
    double sigma = 0.2;
    double t = 0.3846;
    double mu = 0.14;
    int steps = 200;
    double epsilon = 0.01;

    double greeks[5][31];
    double greeks_bs[5][31];

    int s0_increment = 2;
    s0 = 20;
    for(int i = 0; i <= 30; i++){

        double price = s0 + i * s0_increment;
        double curr_price = OptionPricing::callOptionEuropeanBinomial("d", price, k, r, sigma, t, steps);
        double curr_price_bs = OptionPricing::callOptionPriceBS(r,sigma, t, price, k);
        // delta
        greeks[0][i] =  (OptionPricing::callOptionEuropeanBinomial("d", price + epsilon, k, r, sigma, t, steps) - curr_price)/epsilon;
        greeks_bs[0][i] =  (OptionPricing::callOptionPriceBS(r,sigma, t, price + epsilon, k) - curr_price_bs)/epsilon;


        // theta
        greeks[1][i] = (OptionPricing::callOptionEuropeanBinomial("d", price, k, r, sigma, t + epsilon, steps) - curr_price)/epsilon;
        greeks_bs[1][i] =  (OptionPricing::callOptionPriceBS(r,sigma, t + epsilon, price, k) - curr_price_bs)/epsilon;

        // gamma
        greeks[2][i] = (OptionPricing::callOptionEuropeanBinomial("d", price + 1 * 2, k, r, sigma, t, steps)
                - 2 * OptionPricing::callOptionEuropeanBinomial("d", price + 1, k, r, sigma, t, steps) + curr_price)/1;

        greeks_bs[2][i] = ((OptionPricing::callOptionPriceBS(r, sigma , t, price + 1 * 2, k) -
                            2 * OptionPricing::callOptionPriceBS( r, sigma , t, price + 1, k) + curr_price_bs))/(1*1);


        // Vega
        greeks[3][i] = (OptionPricing::callOptionEuropeanBinomial("d", price, k, r, sigma + epsilon * epsilon, t, steps) -
                curr_price) / (epsilon * epsilon);
        greeks_bs[3][i] =  (OptionPricing::callOptionPriceBS(r,sigma + epsilon * epsilon, t, price, k) - curr_price_bs)/(epsilon * epsilon);


        //rho
        greeks[4][i] = (OptionPricing::callOptionEuropeanBinomial("d", price, k, r + epsilon * epsilon, sigma, t, steps) -
                curr_price) / (epsilon * epsilon);
        greeks_bs[4][i] =  (OptionPricing::callOptionPriceBS(r + epsilon * epsilon, sigma, t, price, k) - curr_price_bs)/(epsilon * epsilon);
    }


    Mutils::WriteToCSV2DMatrix(greeks,31,5,"../Data/Qn3.csv");
    Mutils::WriteToCSV2DMatrix(greeks_bs,31,5,"../Data/Qn3_bs.csv");

    s0 = 49;

    int iter = int(t/epsilon);

    double *delta = new double [iter+2];
    double *delta_bs = new double [iter+2];

    t = 0;
    for(int i =0; i<=iter; i++){
        double curr_t = t + i * epsilon;
        double curr_price = OptionPricing::callOptionEuropeanBinomial("d", s0, k, r, sigma, curr_t, steps);
        double curr_price_bs = OptionPricing::callOptionPriceBS(r, sigma, curr_t, s0, k);

        delta[i] =  (OptionPricing::callOptionEuropeanBinomial("d", s0 + epsilon, k, r, sigma, curr_t, steps) - curr_price)/epsilon;
        delta_bs[i] = (OptionPricing::callOptionPriceBS(r, sigma, curr_t, s0 + epsilon, k) - curr_price_bs)/epsilon;
    }
    t = 0.3846;
    double curr_price = OptionPricing::callOptionEuropeanBinomial("d", s0, k, r, sigma, t, steps);
    delta[iter+1] =  (OptionPricing::callOptionEuropeanBinomial("d", s0 + epsilon, k, r, sigma, t, steps) - curr_price)/epsilon;
    Mutils::WriteArrayToCSV(delta,iter+2,"../Data/Qn3_delta_t.csv");

    double curr_price_bs = OptionPricing::callOptionPriceBS(r, sigma, t, s0, k);
    delta_bs[iter+1] = (OptionPricing::callOptionPriceBS(r, sigma, t, s0 + epsilon, k) - curr_price_bs)/epsilon;
    Mutils::WriteArrayToCSV(delta_bs,iter+2,"../Data/Qn3_delta_t_bs.csv");

}

void RunQn4(){
    double r =0.05;
    double sigma = 0.3;
    double k = 100;
    double t = 1;
    double s = 80;
    int steps = 200;

    double * euroPutPrice = new double [11];
    double * americanPutPrice = new double [11];
    double increment = 4;

    for(int i=0;i<=10;i++){
        double stockPrice = s + i * increment;
        euroPutPrice[i] = OptionPricing::putOptionEuropeanBinomial("b", stockPrice, k, r, sigma, t, steps);
        americanPutPrice[i] = OptionPricing::putOptionAmericanBinomial("b", stockPrice, k, r, sigma, t, steps);
    }
    Mutils::WriteArrayToCSV(euroPutPrice,11,"../Data/Qn4EuroPut.csv");
    Mutils::WriteArrayToCSV(americanPutPrice,11,"../Data/Qn4AmericanPut.csv");

}

void RunQn5(){
    double r = 0.05;
    double sigma = 0.24;
    double s0 = 32;
    double k = 30;
    double t = 0.5;

    vector<int> nSize {10, 15, 20, 40, 70, 80, 100, 200, 500};

    cout<<"Option Price Using black sholes: " <<OptionPricing::callOptionPriceBS(r,sigma, t,s0,k)<<endl;

    double *triA = new double [9];
    double *triB = new double [9];

    int counter = 0;

    for(auto &i : nSize) {

        triA[counter] = OptionPricing::callOptionEuropeanTrinomial("a",s0,k, r, sigma, t, i);
        triB[counter] = OptionPricing::callOptionEuropeanTrinomial("b",s0,k, r, sigma, t, i);
        counter++;
    }
    Mutils::WriteArrayToCSV(triA,counter,"../Data/Qn5A.csv");
    Mutils::WriteArrayToCSV(triB,counter,"../Data/Qn5B.csv");

}

void RunQn6(){
    double r = 0.05;
    double sigma = 0.24;
    double s0 = 32;
    double k = 30;
    double t = 0.5;
    double N = 1000;
    int base1 = 5;
    int base2 = 7;

    cout <<"Option Price using Halton's Low Sequence: " << OptionPricing::callOptionEuropeanLDS(s0,k,r,sigma,t,N, base1, base2)<<endl;
    cout <<"Option Price using Black Sholes as Comparision: " <<OptionPricing::callOptionPriceBS(r,sigma,t,s0,k)<<endl;
}

int main(){

    cout << "#################################### Qn 1 ###################################" << endl;
    RunQn1();
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 2 ###################################" << endl;
    RunQn2();
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 3 ###################################" << endl;
    RunQn3();
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 4 ###################################" << endl;
    RunQn4();
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 5 ###################################" << endl;
    RunQn5();
    cout << endl;
    cout << "#################################### Qn 6 ###################################" << endl;
    RunQn6();
    cout << endl;
    return 0;
}

