#include "EquationSolver.h"
#include <limits>

const double eps = std::numeric_limits<double>::epsilon();

class F : public Function{
public:
    double operator () (const double &x) const{
        return x-tan(x);
    }
    double diff (const double &x) const{
        return 1 - 1.0/(cos(x)*cos(x));
    }
} f;

int main(){
    NewtonSolver solver1(f, 4.5, eps, 50);
    double ans1 = solver1.solve();
    NewtonSolver solver2(f, 7.7, eps, 50);
    double ans2 = solver2.solve();
    return 0;
}