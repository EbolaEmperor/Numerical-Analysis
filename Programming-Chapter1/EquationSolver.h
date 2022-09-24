#ifndef _EQUATION_SOLVER_H_
#define _EQUATION_SOLVER_H_

#include <cmath>
#include <iostream>

class Function{
public:
    virtual double operator () (const double &x) const{
        return 0.0;
    }
    virtual double diff (const double &x) const{
        return 0.0;
    }
};

Function baseF;

class EquationSolver{
public:
    virtual double solve(){
        return 0.0;
    }
};

class BisectionSolver : public EquationSolver{
private:
    double a, b, delta, eps;
    Function & f;
    int M;

public:
    BisectionSolver(): f(baseF){}
    BisectionSolver(Function & f, double a, double b, double delta, double eps, int M):
        f(f), a(a), b(b), delta(delta), eps(eps), M(M) {}
    double solve(){
        int k;
        double u = f(a), v = f(b), w, h, c;
        for(k = 1; k <= M; k++){
            h = b - a;
            c = a + h / 2;
            w = f(c);
            if( fabs(h) < delta || fabs(w) < eps )
                break;
            else if( w > 0 && u < 0 || w < 0 && u > 0 )
                b = c, v = w;
            else
                a = c, u = w;
        }
#ifndef SILENCE
        std::cerr << "------------------Bisection Method---------------------" << std::endl;
        std::cerr << "Total Step: " << k << std::endl;
        std::cerr << "Approxiate Root: " << c << std::endl;
        std::cerr << "Function Value: " << f(c) << std::endl;
        std::cerr << "Final Interval Length: " << b-a << std::endl;
#endif
        return c;
    }
};

class NewtonSolver : public EquationSolver{
private:
    double x0, eps;
    Function & f;
    int M;

public:
    NewtonSolver(): f(baseF) {}
    NewtonSolver(Function & f, double x0, double eps, int M):
        f(f), x0(x0), eps(eps), M(M) {}
    double solve(){
        int k;
        double x = x0, u, v;
        for(k = 0; k <= M; k++){
            u = f(x);
            if( fabs(u) < eps ) break;
            v = u / f.diff(x);
            if( x == x - v ) break;
            x = x - v;
        }
#ifndef SILENCE
        std::cerr << "-------------------Newton Method----------------------" << std::endl;
        std::cerr << "Total Step: " << k << std::endl;
        std::cerr << "Initial Point: " << x0 << std::endl;
        std::cerr << "Approxiate Root: " << x << std::endl;
        std::cerr << "Function Value: " << f(x) << std::endl;
        std::cerr << "Derivative Value: " << f.diff(x) << std::endl;
#endif
        return x;
    }
};

class SecandSolver : public EquationSolver{
private:
    double x0, x1, delta, eps;
    Function & f;
    int M;

public:
    SecandSolver(): f(baseF){}
    SecandSolver(Function & f, double x0, double x1, double delta, double eps, int M):
        f(f), x0(x0), x1(x1), delta(delta), eps(eps), M(M) {}
    double solve(){
        int k;
        double u = f(x1), v = f(x0), s;
        for(k = 2; k <= M; k++){
            if( fabs(u) > fabs(v) ){
                std::swap(x0, x1);
                std::swap(u, v);
            }
            s = (x1 - x0) / (u - v);
            x0 = x1;
            v = u;
            x1 = x1 - u * s;
            u = f(x1);
            if( fabs(x1 - x0) < delta || fabs(u) < eps )
                break;
        }
#ifndef SILENCE
        std::cerr << "-------------------Secand Method----------------------" << std::endl;
        std::cerr << "Total Step: " << k << std::endl;
        std::cerr << "Approxiate Root: " << x1 << std::endl;
        std::cerr << "Function Value: " << f(x1) << std::endl;
        std::cerr << "|x(n) - x(n-1)|: " << fabs(x1-x0) << std::endl;
#endif
        return x1;
    }
};

#endif