// 此程序用于测试高精度分数运算下的Hermite插值

#include <iostream>
#include <cmath>
#include "fraction.h"
#include "interpolation.h"

int main(){
    const int n = 5;
    const int xvalues[] = {0, 1, 1, 3, 3};
    const int sp1values[] = {1, 2, -1, 0, 0};

    std::vector<fraction> x(xvalues, xvalues + n);
    std::vector<fraction> sp1(sp1values, sp1values + n);

    T_HermiteInterpolation<fraction> poly1(x, sp1);

    std::cerr << poly1 << std::endl;
    std::cerr << poly1.standardize() << std::endl;
    return 0;
}