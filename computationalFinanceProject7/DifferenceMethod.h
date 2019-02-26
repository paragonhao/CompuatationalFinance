//
// Created by paragonhao on 25/2/19.
//

#ifndef COMPUTATIONALFINANCEPROJECT7_DIFFERENCEMETHOD_H
#define COMPUTATIONALFINANCEPROJECT7_DIFFERENCEMETHOD_H


class DifferenceMethod {

public:
    static double getPUEFD(double dt, double sigma, double delta_X, double r);

    static double getPMEFD(double dt, double sigma, double delta_X, double r);

    static double getPDEFD(double dt, double sigma, double delta_X, double r);

    static void EFDSolver(double currPrice, int deltaFactor);
};


#endif //COMPUTATIONALFINANCEPROJECT7_DIFFERENCEMETHOD_H
