#include <iostream>
#include <cmath>
#include "interpolation.h"

int main(){
    const int n = 10;
    const double xvalues[] = {0, 0, 3, 3, 5, 5, 8, 8, 13, 13};
    const double fvalues[] = {0, 75, 225, 77, 383, 80, 623, 74, 993, 72};

    std::vector<double> x(xvalues, xvalues + n);
    std::vector<double> f(fvalues, fvalues + n);

    HermiteInterpolation hpoly(x, f);
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