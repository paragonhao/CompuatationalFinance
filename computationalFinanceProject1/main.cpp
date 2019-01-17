#include <iostream>
#include <cmath>
#include <cmath>
#include <chrono>
#include "RandomGenerator.h"
#include "Mutils.h"

using namespace std;
using namespace std::chrono;


void RunQn1(int size, long double *arr){
    long double builtRNG[size];
    //  Find mean with the generated random numbers
    double mean = Mutils::Mean(arr, size);

    // Find standard deviation with the generated random numbers
    double sd = Mutils::StDev(arr, size);

    // Genrate the Random Numbers using the built in function of C++
    for(int i =0;i<size;i++){
        builtRNG[i] = ((double) rand() / (RAND_MAX)) ;
    }

    double meanBuiltIn = Mutils::Mean(arr, size);
    double sdBuiltIn = Mutils::StDev(builtRNG, size);

    cout <<"Mean of randomly generated number: " << mean << endl;
    cout <<"Standard Deviation of randomly generated number: " << sd <<endl;
    cout <<"Mean of randomly generated number Using Built in Function: " << meanBuiltIn <<endl;
    cout <<"Standard Deviation of randomly generated number Using Built in Function: " << sdBuiltIn <<endl;
    cout <<"Comment:  LGM random number generator is similar to the random number generator"<<endl;
    cout << endl;
}

void RunQn2(int size, long double *arr){
    long double r[size];

    for(int i=0; i<size; i++){
        if(arr[i] <=0.3 && arr[i] >=0 ){
            r[i] = -1;
        }else if(arr[i]>0.3 && arr[i]<=0.65){
            r[i] = 0;
        }else if(arr[i]>0.65 && arr[i]<=0.85){
            r[i] = 1;
        }else if(arr[i]>0.85 && arr[i]<=1){
            r[i] = 2;
        }
    }

    double meanQn2 = Mutils::Mean(r, size);
    double sdQn2 = Mutils::StDev(r, size);
    cout << "Mean: " <<meanQn2 <<endl;
    cout << "Standard Deviation: "<< sdQn2 <<endl;
    Mutils::WriteToCSV(r, size, "../Data/Q2Data.csv");
    cout << endl;
}

void RunQn3(){
    int binoSize = 1000;
    int n = 44;
    double p = 0.64;

    long double * binoArr = RandomGenerator::rbinom(binoSize, n, p);

    Mutils::WriteToCSV(binoArr, binoSize, "../Data/Q3Data.csv");

    // Compute the probability that P(X>=40)
    int probCounter = 0;
    for(int i=0; i< binoSize; i++){
        probCounter += (binoArr[i] >= 40) ? 1 : 0;
    }

    cout << "Probability that P(X>=40) is : " << probCounter/double(binoSize) << endl;
    cout << "Using R to calculate the exact number: 4.823664e-05. Which is quite similar. " << endl;
}

void RunQn4(int size, long double arr[]){

    long double *expArr = RandomGenerator::rexp(size, arr);

    //Compute P(X>=1) and P(X>=4)
    int x_1=0, x_4=0;
    for(int i=0; i<size; i++){
        x_1 += (expArr[i] >= 1 ) ? 1 : 0;
        x_4 += (expArr[i] >= 4 ) ? 1 : 0;
    }

    Mutils::WriteToCSV(expArr, size, "../Data/Q4Data.csv");

    cout << "Probability that P(X>=1) is : " << x_1/double(size) << endl;
    cout << "Probability that P(X>=4) is : " << x_4/double(size) << endl;

    cout << "Mean of the simulated Exponentail  distribution is: " << Mutils::Mean(expArr, size) << endl;
    cout << "Standard deviation of the simulated Exponentail  distribution is: " << Mutils::StDev(expArr, size) << endl;
}

void RunQn5(int size){

    // generate a common uniform distribution to be used for two algorithms
    long double *arr = RandomGenerator::runif(size);

    // Simulate Standard Normal Distribution Using Box-Muller
    auto start1 = high_resolution_clock::now();

    long double *normArr = RandomGenerator::boxmuller(arr, size);

    auto stop1 = high_resolution_clock::now();
    auto duration1 = duration_cast<microseconds>(stop1 - start1);

    cout << "Simulation Standard Normal Distribution Using Box-Muller: " << endl;
    cout << "Mean: " << Mutils::Mean(normArr, size) << endl;
    cout << "Standard deviation: " << Mutils::StDev(normArr, size) << endl;
    cout << endl;


    // Simulate Standard Normal Distribution Using Polar-Marsaglia
    auto start2 = high_resolution_clock::now();

    long double *normArrPM = RandomGenerator::polarmarsaglia(arr, size);

    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>(stop2 - start2);

    // Number of random variables generated using polar-marsaglia
    int arrSize = int(normArrPM[0]);

    long double *normPM = &normArrPM[1];

    cout << "Simulation Standard Normal Distribution Using Polar-Marsaglia: " << endl;
    cout << "Number of values chosen: " << arrSize << endl;
    cout << "Mean: " << Mutils::Mean(normPM, arrSize) << endl;
    cout << "Standard deviation: " << Mutils::StDev(normPM, arrSize) << endl;
    cout << endl;
    cout << "##################### Comparsion ###################################"<<endl;
    cout << "Simulation with data size: " << size << endl;
    cout << "Time taken for Box-Muller: " << duration1.count() << endl;
    cout << "Time taken for Polar-Marsaglia: " << duration2.count() << endl;
    cout << "Clearly Polar-Marsaglia is faster than Box-Muller" <<endl;
    cout << "####################################################################"<<endl;
    cout << endl;

}


int main() {

    int size = 10000;
    long double * arr = RandomGenerator::runif(size);
    cout << "######################### QN1 ################################"<<endl;
    RunQn1(size, arr);
    cout << "######################### QN2 ################################"<<endl;
    RunQn2(size, arr);
    cout << "######################### QN3 ################################"<<endl;
    RunQn3();
    cout << "######################### QN4 ################################"<<endl;
    RunQn4(size, arr);
    cout << "######################### QN5 ################################"<<endl;
    size = 5000;
    RunQn5(size);

    return 0;
}

