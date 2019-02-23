//
// Created by paragonhao on 22/2/19.
//

#include <cmath>
#include "JDDefaultOption.h"

double JDDefaultOption::getAPR(const double & r0, const double & delta, const double & lamda2){
    return r0 + delta * lamda2;
}

double JDDefaultOption::getPMT(const double & L0, const double & r, const double & T){
    double denominator = 1 - (1 / pow(1 + r, int(T*12)));
    double numerator = L0 * r;

    return numerator / denominator;
}

double JDDefaultOption::getAinLoanFunc(const double & PMT, const double & r){
    return PMT /r;
}

double JDDefaultOption::getBinLoanFunc(const double & PMT, const double & r, const double & T){
    double denominator = r * pow(1 + r, int(T * 12));
    return PMT/denominator;
}