#include <iostream>
#include <cmath>
#include "interpolation.h"

int main(){
    const int n = 7;
    const double xvalues[] = {0.0, 6.0, 10.0, 13.0, 17.0, 20.0, 28.0};
    const double sp1values[] = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    const double sp2values[] = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};

    std::vector<double> x(xvalues, xvalues + n);
    std::vector<double> sp1(sp1values, sp1values + n);
    std::vector<double> sp2(sp2values, sp2values + n);

    NewtonInterpolation poly1(x, sp1);
    NewtonInterpolation poly2(x, sp2);


    NewtonPolynomial::setOutput(NewtonPolynomial::OUTPUT_NORMAL);
    std::cerr << poly1 << std::endl << std::endl;
    std::cerr << poly2 << std::endl << std::endl;
    std::cerr << poly1(43) << " " << poly2(43) << std::endl; 

    // NewtonPolynomial::setOutput(NewtonPolynomial::OUTPUT_LATEX);
    // std::cerr << poly1 << std::endl << std::endl;
    // std::cerr << poly2 << std::endl << std::endl;
    // NewtonPolynomial::setOutput(NewtonPolynomial::OUTPUT_TIKZ);
    // std::cerr << poly1 << std::endl << std::endl;
    // std::cerr << poly2 << std::endl << std::endl;

    return 0;
}