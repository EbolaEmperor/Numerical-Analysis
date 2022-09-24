#include "EquationSolver.h"
#include <limits>

const double pi = acos(-1);
const double eps = 1e-6;

const double L = 10, r = 1, V = 12.4;

class F : public Function{
public:
    double operator () (const double &h) const{
        return L * ( 0.5 * pi * r * r - r * r * asin(h/r) - h * sqrt(r * r - h * h ) ) - V;
    }
    double diff (const double &h) const{
        return L * ( - r * 1.0 / sqrt(1.0 - (h / r) * (h / r)) - sqrt(r * r - h * h ) + 2 * h * h / sqrt(r * r - h * h) );
    }
} f;

int main(){
    NewtonSolver sol1(f, 0, eps, 100);
    double ans1 = sol1.solve();
    BisectionSolver sol2(f, 0, 1, eps, eps, 100);
    double ans2 = sol2.solve();
    SecandSolver sol3(f, 0, 0.5, eps, eps, 100);
    double ans3 = sol3.solve();
    return 0;
}