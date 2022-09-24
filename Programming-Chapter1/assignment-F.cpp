#include "EquationSolver.h"
#include <limits>

const double pi = acos(-1);
const double eps = std::numeric_limits<double>::epsilon();

class F : public Function{
private:
    double A, B, C, E;
public:
    F(){}
    F(double l, double h, double D, double b){
        double sinb = sin(b/180*pi);
        double cosb = cos(b/180*pi);
        double tanb = tan(b/180*pi);
        A = l*sinb;
        B = l*cosb;
        C = (h+0.5*D)*sinb - 0.5*D*tanb;
        E = (h+0.5*D)*cosb-0.5*D;
    }
    double operator () (const double &x) const{
        double sinx = sin(x);
        double cosx = cos(x);
        return A*sinx*cosx + B*sinx*sinx - C*cosx - E*sinx;
    }
    double diff(const double &x) const{
        double sinx = sin(x);
        double cosx = cos(x);
        return A*(-sinx*sinx+cosx*cosx) + 2*B*cosx*sinx + C*sinx - E*cosx;
    }
};

void test1(){
    F f(89, 49, 55, 11.5);
    NewtonSolver sol(f, 33.0/180.0*pi, eps, 100);
    double ans = sol.solve();
    std::cerr << "ans = " << ans/pi*180 << " degree" << std::endl << std::endl;
}

void test2(){
    F f(89, 49, 30, 11.5);
    NewtonSolver sol(f, 33.0/180.0*pi, eps, 100);
    double ans = sol.solve();
    std::cerr << "ans = " << ans/pi*180 << " degree" << std::endl << std::endl;
}

void test3(){
    F f(89, 49, 55, 11.5);
    SecandSolver sol(f, 0, 5.0/180.0*pi, eps, eps, 100);
    double ans = sol.solve();
    std::cerr << "ans = " << ans/pi*180 << " degree" << std::endl << std::endl;
}

int main(){
    test1();
    test2();
    test3();
    return 0;
}