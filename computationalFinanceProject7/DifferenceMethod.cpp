//
// Created by paragonhao on 25/2/19.
//

#include "DifferenceMethod.h"

double DifferenceMethod::getPUEFD(double dt, double sigma, double delta_X, double r){
    return dt * ((sigma * sigma)/ (2 * delta_X * delta_X) + (r - 0.5 * sigma * sigma)/ (2 * delta_X));
}

double DifferenceMethod::getPMEFD(double dt, double sigma, double delta_X, double r){
    return 1 - (dt * sigma * sigma / (delta_X * delta_X)) - r * dt;
}

double DifferenceMethod::getPDEFD(double dt, double sigma, double delta_X, double r){
    return dt * ((sigma * sigma)/ (2 * delta_X * delta_X) - (r - 0.5 * sigma * sigma)/ (2 * delta_X));
}