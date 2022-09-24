#include "EquationSolver.h"
#include <limits>

const double pi = acos(-1);
const double eps = std::numeric_limits<double>::epsilon();

//--------------------test1--------------------------------

class Func1 : public Function{
public:
    double operator () (const double &x) const{
        return sin(x/2) - 1.0;
    }
} func1;

void test1(){
    SecandSolver solver(func1, 0, pi/2, eps, eps, 100);
    double ans = solver.solve();
}

//--------------------test2--------------------------------

class Func2 : public Function{
public:
    double operator () (const double &x) const{
        return exp(x) - tan(x);
    }
} func2;

void test2(){
    SecandSolver solver(func2, 1.0, 1.4, eps, eps, 100);
    double ans = solver.solve();
}

//--------------------test3--------------------------------

class Func3 : public Function{
public:
    double operator () (const double &x) const{
        return x*x*x - 12*x*x + 3*x + 1;
    }
} func3;

void test3(){
    SecandSolver solver(func3, 0, -0.5, eps, eps, 100);
    double ans = solver.solve();
}

int main(){
    std::cerr << "sin(x/2) - 1.0" << std::endl;
    test1();
    std::cerr << "exp(x) - tan(x)" << std::endl;
    test2();
    std::cerr << "x*x*x - 12*x*x + 3*x + 1" << std::endl;
    test3();
    return 0;
}