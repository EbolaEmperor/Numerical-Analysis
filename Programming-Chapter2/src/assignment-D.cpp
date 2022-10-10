#include <iostream>
#include <cmath>
#include "interpolation.h"

int main(){
    const int n = 5;
    const double xvalues[] = {0, 3, 5, 8, 13};
    const double fvalues[] = {0, 225, 383, 623, 993};
    const double dfvalues[] = {75, 77, 80, 74, 72};

    std::vector<double> x(xvalues, xvalues + n);
    std::vector<double> f(fvalues, fvalues + n);
    std::vector<double> df(dfvalues, dfvalues + n);

    HermiteInterpolation hpoly(x, f, df);
    Polynomial poly = hpoly.standardize();
    Polynomial dpoly = poly.diff();

    Polynomial::setOutput(Polynomial::OUTPUT_LATEX);
    std::cout << poly << std::endl << std::endl;
    std::cout << poly(10) << std::endl << std::endl;
    std::cout << dpoly << std::endl << std::endl;

    // 若需要输出 Tikz可识别 的格式， 请将下面的注释取消
    // Polynomial::setOutput(Polynomial::OUTPUT_TIKZ);
    // std::cout << poly << std::endl << std::endl;
    // std::cout << poly(10) << std::endl << std::endl;
    // std::cout << dpoly << std::endl << std::endl;

    return 0;
}