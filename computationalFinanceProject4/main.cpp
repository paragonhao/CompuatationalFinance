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

    vector<int> nSize {10, 20, 40, 80, 100, 200};
    vector<string> method {"a", "b", "c", "d"};
    cout << "Option price using Black Sholes: " << OptionPricing::callOptionPriceBS(r,sigma,t, s0, k) << endl;

    for(int j=0; j<4; j++){
        for(int i=0; i< 6; i++){
            cout << "Option Price for method ("<< method[j]<< "), n = "<< nSize[i] <<" is: " << OptionPricing::callOptionEuropeanBinomial(method[j], s0, k, r, sigma, t, nSize[i]) << endl;
        }
        cout << "####################################################################################" << endl;
    }
}

void RunQn2(){
    double rf = 0.02;
    double price = 1141.42;
    double k = (int(price)/10) * 10.0 * 1.1;

    ifstream  data("../Data/GOOG_2013_01_2018_01.csv");


    string line;
    int counter = 0;
    long double *adjustedPrice = new long double[1500];

    while(getline(data,line))
    {
        long double temp = stold(line);
        adjustedPrice[counter] = temp;
        counter ++;
    }

    vector<int> yearEndpoints = {251, 502, 753, 1006};
    long double ret[4];

    int prev = 0;
    for(int i=0; i<6; i++){
        if( i< 4){
            int pos = yearEndpoints[i];
            ret[i] = (adjustedPrice[pos] - adjustedPrice[prev])/adjustedPrice[prev];
            prev = pos;
        }
    }
    long double sigma = Mutils::StDev(ret, 4);
    cout << "Price is 58.8 on Yahoo Finance"<< endl;
    cout << "Estimate the option price on Jan 2020: " <<OptionPricing::callOptionEuropeanBinomial("d", price, k, rf, sigma, 1, 200)<< endl;
    cout << "Estimate the option price with Volatility equal to 20.4%: " <<OptionPricing::callOptionEuropeanBinomial("d", price, k, rf, 0.204, 1, 200)<< endl;

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

        double price = s0 + i * 2;
        double curr_price = OptionPricing::callOptionEuropeanBinomial("d", price, k, r, sigma, t, steps);
        double curr_price_bs = OptionPricing::callOptionPriceBS(r,sigma, t, s0, k);
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
    cout << "output data"<<endl;

    Mutils::WriteToCSV2DMatrix(greeks,31,5,"../Data/Qn3.csv");




}

int main() {

    cout << "#################################### Qn 1 ###################################" << endl;
    //RunQn1();
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 2 ###################################" << endl;
    //RunQn2();
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 3 ###################################" << endl;
    RunQn3();
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 4 ###################################" << endl;
//    RunQn4(seed);
    cout << "#############################################################################" << endl;
    cout << endl;
    cout << "#################################### Qn 5 ###################################" << endl;
//    RunQn5(seed);
    cout << endl;
    return 0;
}

