#include <iostream>
#include <cmath>
#include "interpolation.h"

class Func : public Function{
public:
    double operator () (const double &x) const{
        return 1.0 / (1.0 + x * x);
    }
} func;

int main(){
    std::vector<NewtonInterpolation> poly;

    for(int n = 2; n <= 8; n += 2){
        std::vector<double> x;
        for(int i = 0; i <= n; i++)
            x.push_back(-5.0 + 10.0 * i / n);
        poly.push_back(NewtonInterpolation(func, x));
    }

    NewtonInterpolation::setOutput(NewtonInterpolation::OUTPUT_NORMAL);
    for(auto p : poly) std::cout << p << std::endl << std::endl;

    NewtonInterpolation::setOutput(NewtonInterpolation::OUTPUT_LATEX);
    for(auto p : poly) std::cout << p << std::endl << std::endl;

    NewtonInterpolation::setOutput(NewtonInterpolation::OUTPUT_TIKZ);
    for(auto p : poly) std::cout << p << std::endl << std::endl;

    return 0;
}