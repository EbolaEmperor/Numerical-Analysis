#include <iostream>
#include <cmath>
#include "interpolation.h"

class Func : public Function{
public:
    double operator () (const double &x) const{
        return 6*pow(x,6) - pow(x,4) + 5*pow(x,3) - 3*pow(x,2) + x - 9;
    }
} func;

int main(){
    std::vector<double> x;
    x.push_back(1);
    x.push_back(2);
    x.push_back(3);
    x.push_back(4);
    x.push_back(5);
    x.push_back(6);
    x.push_back(7);

    NewtonInterpolation p(func, x);
    for(int i = 0; i < 10; i++){
        double vx = (double)rand()/RAND_MAX*10-5;
        std::cout << vx << " " << func(vx) - p(vx) << std::endl;
    }
    return 0;
}