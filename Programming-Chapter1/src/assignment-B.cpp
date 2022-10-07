#include "EquationSolver.h"
#include <limits>

const double pi = acos(-1);
const double eps = std::numeric_limits<double>::epsilon();

//--------------------test1--------------------------------

class Func1 : public Function{
public:
    double operator () (const double &x) const{
        return 1.0/x-tan(x);
    }
} func1;

void test1(){
    BisectionSolver solver(func1, 0, pi/2, eps, eps, 100);
    double ans = solver.solve();
}

//--------------------test2--------------------------------

class Func2 : public Function{
public:
    double operator () (const double &x) const{
        return 1.0/x-pow(2.0,x);
    }
} func2;

void test2(){
    BisectionSolver solver(func2, 0, 1, eps, eps, 100);
    double ans = solver.solve();
}

//--------------------test3--------------------------------

class Func3 : public Function{
public:
    double operator () (const double &x) const{
        return pow(2.0, -x) + exp(x) + 2.0*cos(x) - 6.0;
    }
} func3;

void test3(){
    BisectionSolver solver(func3, 1, 3, eps, eps, 100);
    double ans = solver.solve();
}

//--------------------test4--------------------------------

class Func4 : public Function{
public:
    double operator () (const double &x) const{
        return (x*x*x + 4*x*x + 3*x + 5) / (2*x*x - 9*x*x + 18*x - 2);
    }
} func4;

void test4(){
    BisectionSolver solver(func4, 0, 4, eps, eps, 100);
    double ans = solver.solve();
}

int main(){
    std::cerr << "1.0/x-tan(x)" << std::endl;
    test1();
    std::cerr << "1.0/x-pow(2.0,x)" << std::endl;
    test2();
    std::cerr << "pow(2.0, -x) + exp(x) + 2.0*cos(x) - 6.0" << std::endl;
    test3();
    std::cerr << "(x*x*x + 4*x*x + 3*x + 5) / (2*x*x - 9*x*x + 18*x - 2)" << std::endl;
    test4();
    return 0;
}