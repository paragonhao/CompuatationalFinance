//
// Created by paragonhao on 25/2/19.
//

#ifndef COMPUTATIONALFINANCEPROJECT7_DIFFERENCEMETHOD_H
#define COMPUTATIONALFINANCEPROJECT7_DIFFERENCEMETHOD_H

#include <string>

using namespace std;

class DifferenceMethod {

public:
    static double getPUEFD(double dt, double sigma, double delta_X, double r);

    static double getPMEFD(double dt, double sigma, double delta_X, double r);

    static double getPDEFD(double dt, double sigma, double delta_X, double r);

    static double getPDIFD(double dt, double sigma, double delta_X, double r);

    static double getPMIFD(double dt, double sigma, double delta_X, double r);

    static double getPUIFD(double dt, double sigma, double delta_X, double r);

    static double getPUCNFD(double dt, double sigma, double delta_X, double r);

    static double getPMCNFD(double dt, double sigma, double delta_X, double r);

    static double getPDCNFD(double dt, double sigma, double delta_X, double r);

    static void EFDEuroPutSolver(double currPrice, int deltaFactor);

    static void IFDEuroPutSolver(double currPrice, int deltaFactor);

    static void CNFDEuroPutSolver(double currPrice, int deltaFactor);

    static void GeneralisationOptionPriceSolver(double currPrice, double ds, string method, string type);

};


#endif //COMPUTATIONALFINANCEPROJECT7_DIFFERENCEMETHOD_H
