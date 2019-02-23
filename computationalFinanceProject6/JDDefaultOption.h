//
// Created by paragonhao on 22/2/19.
//

#ifndef COMPUTATIONALFINANCEPROJECT6_JDDEFAULTOPTION_H
#define COMPUTATIONALFINANCEPROJECT6_JDDEFAULTOPTION_H


class JDDefaultOption {

public:
    static double getAPR(const double & r0, const double & delta, const double & lamda2);

    static double getAinLoanFunc(const double &PMT, const double &r);

    static double getBinLoanFunc(const double &PMT, const double &r, const double &T);

    static double getPMT(const double &L0, const double &r, const double &T);
};


#endif //COMPUTATIONALFINANCEPROJECT6_JDDEFAULTOPTION_H
