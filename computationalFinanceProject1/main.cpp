#include <iostream>
#include <math.h>
#include <cmath>
#include <chrono>
#include "Utils.h"
#include "RandomGenerator.h"

using namespace std;
using namespace std::chrono;


void RunQn1(int size, long double *arr){
    long double builtRNG[size];
    //  Find mean with the generated random numbers
    double mean = Mean(arr, size);

    // Find standard deviation with the generated random numbers
    double sd = StDev(arr, size);

    // Genrate the Random Numbers using the built in function of C++
    for(int i =0;i<size;i++){
        builtRNG[i] = ((double) rand() / (RAND_MAX)) ;
    }

    double meanBuiltIn = Mean(arr, size);
    double sdBuiltIn = StDev(builtRNG, size);

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

    double meanQn2 = Mean(r, size);
    double sdQn2 = StDev(r, size);
    cout << "Mean: " <<meanQn2 <<endl;
    cout << "Standard Deviation: "<< sdQn2 <<endl;
    WriteToCSV(r, size, "../Data/Q2Data.csv");
    cout << endl;
}

void RunQn3(){
    int binoSize = 1000;
    int n = 44;
    double p = 0.64;

    long double * binoArr = RandomGenerator::rbinom(binoSize, n, p);

    WriteToCSV(binoArr, binoSize, "../Data/Q3Data.csv");

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

    WriteToCSV(expArr, size, "../Data/Q4Data.csv");

    cout << "Probability that P(X>=1) is : " << x_1/double(size) << endl;
    cout << "Probability that P(X>=4) is : " << x_4/double(size) << endl;

    cout << "Mean of the simulated Exponentail  distribution is: " << Mean(expArr, size) << endl;
    cout << "Standard deviation of the simulated Exponentail  distribution is: " << StDev(expArr, size) << endl;
}

void RunQn5(int size){
    long double *arr = RandomGenerator::runif(size);

    double normArr[size];

    // Simulate Standard Normal Distribution Using Box-Muller
    auto start1 = high_resolution_clock::now();

    for(int i=0; i< size; i+=2){
        // calculate z1
        normArr[i] = sqrt(-2 * log(arr[i])) * cos(2 * M_PI * arr[i+1]);

        // calculate z2
        normArr[i+1] = sqrt(-2 * log(arr[i])) * sin(2 * M_PI * arr[i+1]);
    }
    auto stop1 = high_resolution_clock::now();
    auto duration1 = duration_cast<microseconds>(stop1 - start1);

    cout << "Simulation Standard Normal Distribution Using Box-Muller: " << endl;
    cout << "Mean: " << Mean(normArr, size) << endl;
    cout << "Standard deviation: " << StDev(normArr, size) << endl;
    cout << endl;


    // Simulate Standard Normal Distribution Using Polar-Marsaglia
    double normArrPM[size];
    int arrSize = 0;
    double v1 = 0;
    double v2 = 0;
    double w = 0;

    // Verify arr is a uniform distribution
    //cout << Mean(arr, size)<< endl;
    //cout << StDev(arr, size)<< endl;

    auto start2 = high_resolution_clock::now();

    for(int i=0; i<size; i+=2){
        v1 = 2 * arr[i] - 1;
        v2 = 2 * arr[i+1] -1;
        w = v1 * v1  + v2 * v2; // using pow() would significantly increase the execution time.

        if(w <= 1.0 ){
            normArrPM[arrSize] = v1 * sqrt(-2 * log(w) / w);
            normArrPM[arrSize+1] = v2 * sqrt(-2 * log(w) / w);
            arrSize +=2;
        }
    }
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>(stop2 - start2);

    cout << "Simulation Standard Normal Distribution Using Polar-Marsaglia: " << endl;
    cout << "Number of values chosen: " << arrSize << endl;
    cout << "Mean: " << Mean(normArrPM, arrSize) << endl;
    cout << "Standard deviation: " << StDev(normArrPM, arrSize) << endl;
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
    //    Findings for 1(c):
    //    Mean of randomly generated number: 0.501751
    //    Standard Deviation of randomly generated number: 0.290391
    //    Mean of randomly generated number Using Built in Function: 0.501751
    //    Standard Deviation of randomly generated number Using Built in Function: 0.28921
    //
    //    Comment:  LGM random number generator is similar to the random number generator
    long double m = pow(2, 31) -1;

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

