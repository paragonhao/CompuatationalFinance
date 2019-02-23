#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"
#include "JDDefaultOption.h"
#include <random>

double * GenerateQtFunction(int T, double t, int steps){
    double alpha = 0.7;
    double epsilon = 0.95;
    double beta = (epsilon - alpha) / T;

    auto * qt = new double [steps +1];
    qt[0] = alpha + beta * t * 0;
    for(int i=1; i <= steps;i++){
        qt[i] = alpha + beta * t * i;
    }

    return qt;

}

double *  GenerateLFunction(double T, double t, int steps, double lamdaTwo){
    double L0 = 22000;
    double r0 = 0.02;
    double delta = 0.25;


    double R = r0 + delta * lamdaTwo;
    double r = R/12;
    double pmt = JDDefaultOption::getPMT(L0, r, T);
    double a = JDDefaultOption::getAinLoanFunc(pmt, r);
    double b = JDDefaultOption::getBinLoanFunc(pmt,r,T);
    double c = 1 + r;

    auto * LoanArray = new double [steps +1];
    LoanArray[0] = 22000;

    for(int i =1; i <= steps;i ++){
        LoanArray[i] = a - b * pow(c , 12 * t * i);
    }
    LoanArray[steps] = 0;

    return LoanArray;
}

double * JumpDiffusionModelCollatoral(const double & t, const int &steps, int seed, double lamdaOne){
    // Run simulation for one path
    double V0 = 20000;
    double gamma = -0.4;
    double mu = -0.1;
    double sigma = 0.2;

    default_random_engine generator(seed);
    exponential_distribution<double> distribution(lamdaOne);
    // first jump
    double exp1 = distribution(generator);

    auto * value = new double [steps+1];
    auto * dwt = RandomGenerator::wienerProcess(t,steps, seed);
    value[0] = V0;

    for(int i =0; i<steps; i++){
        double currTime = t * i;
        if(currTime > exp1){

            exp1 += distribution(generator);
            value[i] = value[i] * (1 + gamma);
        }
        value[i + 1] = value[i] + value[i] * mu * t + sigma * value[i] * dwt[i];
        if(value[i + 1] < 0.01){
            value[i + 1] =0;
        }
    }
    return value;
}

void RunQ1(){
    double r =0.03;
    double T = 1;
    double s0 = 98;
    double strike = 100;
    double sigma = 0.12;
    double increment = 0.04;

    int num = 1000;
    int steps = 100;
    int sigmaSteps =10;
    double delta = T/steps;

    cout <<"Call Price is: "<<endl;
    for(int i=0; i< sigmaSteps; i++){
        cout << OptionPricing::FixedStrikeLookBackCall(r, T, s0, strike, sigma, num, steps, delta) << endl;
        sigma += increment;
    }
    cout <<"Put Price is: "<<endl;
    sigma = 0.12;
    for(int i=0; i< sigmaSteps; i++){
        cout << OptionPricing::FixedStrikeLookBackPut(r, T, s0, strike, sigma, num, steps, delta) << endl;
        sigma += increment;
    }

}

double * Proj6_2function(const double & lamdaOne, const double & lamdaTwo, const double & T){
    //base parameters;
    double epsilon = 0.95;
    double r0 = 0.02;
    // end of base parameters

    // parameters used in a loop
    int  steps = 1000; // simulate 1000 paths for the process, including the starting  point
    double t = T/steps;
    int nSims = 1000;
    double LQ = 0, LS = 0, VQ = 0, VS = 0;
    double Q = 0, S = 0;
    int defaultCounter = 0;
    double expectedTau = 0;
    double tau = 0;
    int seed =1234;
    // These two values are fixed
    double * LoanArray = GenerateLFunction(T, t, steps, lamdaTwo); // resulting array is steps + 1
    double * Qt = GenerateQtFunction(T, t, steps);   // resulting array is steps + 1
    double totalPayoff = 0;
    //

    // generate 100 simulation path
    for(int k =0; k< nSims; k++){ //############ outer loop ####################################
        // initialize a seed from each iteration
        seed = seed + 500 * k;
        tau = T;
        Q = T;

        default_random_engine generator(seed);
        exponential_distribution<double> distribution(lamdaTwo);

        // generate S, it is the extreme event
        S = distribution(generator);

        // simulate a jump model for the value path
        double * value = JumpDiffusionModelCollatoral(t, steps, seed, lamdaOne);


        for(int i =0; i<= steps; i++){ //############ inner loop ####################################
            double currTime = t * i;
            double QtLt = Qt[i] * LoanArray[i];

            if(currTime > S){

                VS = value[i];
                LS = LoanArray[i];
                tau = S;
                break;
            }else if(value[i] <= QtLt ){
                Q = currTime;
                VQ = value[i];
                LQ = LoanArray[i];
                tau = currTime;
            }
        }

        if(tau < T){
            defaultCounter++;
            expectedTau += tau;

            // Q is exercised first
            if(Q <= S){
                // pay off need to be discounted back from time Q
                totalPayoff += Mutils::max(LQ - epsilon * VQ, 0 ) * exp(-r0 * Q);
            }else{
                totalPayoff += abs(LS - epsilon * VS) * exp(-r0 * S);
            }
        }
    }

    auto * result = new double[3]; // first is Option price, second is probability, third is Expected Tau
    result[0] = totalPayoff/nSims;
    result[1] = double(defaultCounter)/nSims;
    result[2] = expectedTau/defaultCounter;
    return result;
}



int main() {
    cout << "######################################## Qn1 ##################################################" <<endl;
    RunQ1();
    cout << "###############################################################################################" <<endl;
    cout << "######################################## Qn2 ##################################################" <<endl;

    double lamdaTwo = 0.4;
    double lamdaOneIncrement = 0.05;
    double lamdaOne = 0.05;
    cout << "lamda Two is 0.4" << endl;
    cout << "Year lamdaOne Option_Price Probability Expected_Tau"<< endl;

    for(int t = 3; t <= 8; t++){
        for(int i = 0; i<=8;i++){
            double * result = Proj6_2function((lamdaOne + i*lamdaOneIncrement), lamdaTwo, t);
            cout << t<<" "<<(lamdaOne + i*lamdaOneIncrement)<< " "<<result[0]<<" "<< result[1] <<" "<<result[2]<<endl;
        }
    }

    lamdaOne = 0.2;
    double lamdaTwoIncrement = 0.1;
    cout << "lamda One is 0.2" << endl;
    cout << "Year lamdaTwo Option_Price Probability Expected_Tau"<< endl;
    for(int t = 3; t <= 8; t++){
        for(int i = 0; i<=8;i++){
            double * result = Proj6_2function(lamdaOne, i * lamdaTwoIncrement, t);
            cout << t<<" "<<(i * lamdaTwoIncrement)<< " "<<result[0]<<" "<< result[1] <<" "<<result[2]<<endl;
        }
    }



    cout << "###############################################################################################" <<endl;
}