#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include "RandomGenerator.h"
#include "Mutils.h"
#include "OptionPricing.h"
#include "JDDefaultOption.h"
#include "DifferenceMethod.h"
#include <random>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main() {
    cout << "######################################## Qn1 ##################################################" <<endl;
    cout << "Explicit Finite-Difference method: sigma * sqrt(dx)" <<endl;
    cout << "Price, Black sholes, Put Pay off" <<endl;
    int deltaFactor = 1;
    DifferenceMethod::EFDSolver(10, deltaFactor);
    cout << "sigma * sqrt(3 * dx)" <<endl;
    deltaFactor = 3;
    DifferenceMethod::EFDSolver(10, deltaFactor);
    cout << "sigma * sqrt(4 * dx)" <<endl;
    deltaFactor = 4;
    DifferenceMethod::EFDSolver(10, deltaFactor);
    cout << "###############################################################################################" <<endl;



    return 0;
}