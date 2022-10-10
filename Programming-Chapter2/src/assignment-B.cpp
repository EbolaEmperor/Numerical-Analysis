#include <iostream>
#include <cmath>
#include "interpolation.h"

class Func{
public:
    double operator () (const double &x) const{
        return 1.0 / (1.0 + x * x);
    }
} func;

int main(){
    std::vector<NewtonInterpolation> poly;

    for(int n = 2; n <= 8; n += 2){
        std::vector<double> x;
        std::vector<double> f;
        for(int i = 0; i <= n; i++){
            x.push_back(-5.0 + 10.0 * i / n);
            f.push_back(func(x[i]));
        }
        poly.push_back(NewtonInterpolation(x, f));
    }

    NewtonPolynomial::setOutput(NewtonPolynomial::OUTPUT_NORMAL);
    for(auto p : poly) std::cout << p << std::endl << std::endl;

    // 若需要输出 Latex公式 或 Tikz可识别 的格式， 请将下面的注释取消

    // NewtonPolynomial::setOutput(NewtonPolynomial::OUTPUT_LATEX);
    // for(auto p : poly) std::cout << p << std::endl << std::endl;
    // NewtonPolynomial::setOutput(NewtonPolynomial::OUTPUT_TIKZ);
    // for(auto p : poly) std::cout << p << std::endl << std::endl;

    return 0;
}