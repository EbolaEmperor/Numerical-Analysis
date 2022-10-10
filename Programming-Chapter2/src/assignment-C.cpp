#include <iostream>
#include <cmath>
#include "interpolation.h"

const double pi = acos(-1);

class Func{
public:
    double operator () (const double &x) const{
        return 1.0 / (1.0 + 25.0 * x * x);
    }
} func;

int main(){
    std::vector<NewtonInterpolation> poly;

    for(int n = 5; n <= 20; n += 5){
        std::vector<double> x;
        std::vector<double> f;
        for(int i = 1; i <= n; i++){
            x.push_back( cos((2.0 * i - 1.0) / (2.0 * n) * pi) );
            f.push_back(func(x[i-1]));
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