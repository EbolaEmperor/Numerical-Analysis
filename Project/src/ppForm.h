#ifndef _PPFORM_H_
#define _PPFORM_H_

#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <limits>
#include "polynomial.h"
#include "matrix.h"
#include "Function.h"

class ppForm_base{
protected:
    int n;
    std::vector<double> knots;
    std::vector<Polynomial> poly;

public:
    ppForm_base() {}

    ppForm_base(const ppForm_base &rhs):
        n(rhs.n), knots(rhs.knots), poly(rhs.poly) {}
    
    double operator () (const double &x) const{
        if(n == 0)
            error("ppForm operator () :: Visited an uninitialized object.");
        if(x < knots.front() || x > knots.back())
            error("ppForm operator () :: point x are not in the interpolation interval.");
        for(int i = 0; i < n-1; i++){
            if(x <= knots[i+1]) return poly[i](x);
        }
        return poly[n-1](x);
    }
};

class ppForm_linear : public ppForm_base{
public:
    ppForm_linear(const std::vector<double> &x, const std::vector<double> &f){
        if(x.size() != f.size() || x.size() < 2)
            error("[Error] ppForm_linear initializing :: The knots and values are incorrect.");
        n = x.size() - 1;
        for(int i = 0; i < n; i++)
            if(x[i] >= x[i+1])
                error("ppForm_linear initializing :: The knots are not increasing.");
        knots = x;
        for(int i = 0; i < n; i++){
            poly.push_back( Polynomial( (x[i+1]*f[i]-x[i]*f[i+1])/(x[i+1]-x[i]), (f[i+1]-f[i])/(x[i+1]-x[i]) ) );
        }
    }

    ppForm_linear(const std::vector<double> &x, Function & func){
        std::vector<double> f;
        for(auto & xv : x) f.push_back(func(xv));
        (*this) = ppForm_linear(x, f);
    }

    ppForm_linear(const int &n, const double &l, const double &r, Function &func){
        std::vector<double> t(n);
        for(int i = 0; i < n-1; i++)
            t[i] = l + (r-l)*i/(n-1);
        t.back() = r;
        (*this) = ppForm_linear(t, func);
    }
};

class ppForm_cubic : public ppForm_base{
private:
    void basic_init(Matrix &A, ColVector &b, const std::vector<double> &x, const std::vector<double> &f){
        // 端点处函数值条件
        for(int i = 0; i < n; i++){
            for(int j = 0; j < 4; j++){
                A[2*i][4*i+j] = pow(x[i], j);
                A[2*i+1][4*i+j] = pow(x[i+1], j);
            }
            b[2*i] = f[i];
            b[2*i+1] = f[i+1];
        }
        // 端点处一阶、二阶导函数连续条件
        for(int i = 0; i < n-1; i++){
            for(int j = 1; j < 4; j++){
                A[2*n + 2*i][4*i+j] = j * pow(x[i+1], j-1);
                A[2*n + 2*i][4*(i+1)+j] = - j * pow(x[i+1], j-1);
            }
            for(int j = 2; j < 4; j++){
                A[2*n + 2*i+1][4*i+j] = j*(j-1) * pow(x[i+1], j-2);
                A[2*n + 2*i+1][4*(i+1)+j] = - j*(j-1) * pow(x[i+1], j-2);
            }
            b[2*n + 2*i] = 0;
            b[2*n + 2*i+1] = 0;
        }
    }

public:
    ppForm_cubic(const std::vector<double> &x, const std::vector<double> &f, const std::string &bondary){
        for(int i = 0; i < x.size()-1; i++)
            if(x[i] >= x[i+1]) error("ppForm_cubic :: the knots must be strictly increasing!");
        if(f.size() < x.size())
            error("ppForm_cubic :: Too small size of vector<double> f!");
        
        n = x.size() - 1;
        Matrix A(4*n, 4*n);
        ColVector b(4*n);
        basic_init(A, b, x, f);
        if(bondary == "natural"){
            for(int j = 2; j < 4; j++){
                A[4*n-2][j] = j*(j-1) * pow(x[0], j-2);
                A[4*n-1][4*(n-1)+j] = j*(j-1) * pow(x[n], j-2);
            }
            b[4*n-2] = b[4*n-1] = 0;
        } else if(bondary == "complete"){
            if(f.size() < x.size()+2)
                error("ppForm_cubic :: Cannot read derivatives at ends in vector<double> f!");
            for(int j = 1; j < 4; j++){
                A[4*n-2][j] = j * pow(x[0], j-1);
                A[4*n-1][4*(n-1)+j] = j * pow(x[n], j-1);
            }
            b[4*n-2] = f[n+1];
            b[4*n-1] = f[n+2];
        } else if(bondary == "not-a-knot"){
            A[4*n-2][3] = 1;
            A[4*n-2][7] = -1;
            A[4*n-1][4*(n-2)+3] = 1;
            A[4*n-2][4*(n-1)+3] = -1;
            b[4*n-2] = b[4*n-1] = 0;
        } else if(bondary == "second-derivatives-at-end"){
            if(f.size() < x.size()+2)
                error("ppForm_cubic :: Cannot read second-derivatives at ends in vector<double> f!");
            for(int j = 2; j < 4; j++){
                A[4*n-2][j] = j*(j-1) * pow(x[0], j-2);
                A[4*n-1][4*(n-1)+j] = j*(j-1) * pow(x[n], j-2);
            }
            b[4*n-2] = f[n+1];
            b[4*n-1] = f[n+2];
        } else if(bondary == "periodic"){
            if(f.front() != f.back())
                error("ppForm_cubic :: The left value and right value are not the same in boundary 'periodic'!");
            for(int j = 1; j < 4; j++){
                A[4*n-2][j] = j * pow(x[0], j-1);
                A[4*n-2][4*(n-1)+j] = - j * pow(x[n], j-1);
            }
            for(int j = 2; j < 4; j++){
                A[4*n-1][j] = j*(j-1) * pow(x[0], j-2);
                A[4*n-1][4*(n-1)+j] = - j*(j-1) * pow(x[n], j-2);
            }
            b[4*n-2] = b[4*n-1] = 0;
        }
        ColVector coef = solve(A, b);
        std::vector<double> pcoef(4);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < 4; j++)
                pcoef[j] = coef[4*i+j];
            poly.push_back(Polynomial(pcoef));
        }
        knots = x;
    }

    ppForm_cubic(const std::vector<double> &x, const std::vector<double> &f):
        ppForm_cubic(x, f, "natural") {}

    ppForm_cubic(const std::vector<double> &x, Function &f, const std::string &bondary){
        std::vector<double> fv;
        for(auto v : x) fv.push_back(f(v));
        if(bondary == "complete"){
            fv.push_back(f.diff(x.front()));
            fv.push_back(f.diff(x.back()));
        } else if(bondary == "second-derivatives-at-end"){
            fv.push_back(f.diff2(x.front()));
            fv.push_back(f.diff2(x.back()));
        } else if(bondary == "periodic"){
            fv.back() = f(x[0]);
        }
        (*this) = ppForm_cubic(x, fv, bondary);
    }

    ppForm_cubic(const std::vector<double> &x, Function &f):
        ppForm_cubic(x, f, "natural") {}

    ppForm_cubic(const int &n, const double &l, const double &r, Function & func, const std::string &bondary){
        std::vector<double> t(n);
        for(int i = 0; i < n-1; i++)
            t[i] = l + (r-l)*i/(n-1);
        t.back() = r;
        (*this) = ppForm_cubic(t, func, bondary);
    }

    ppForm_cubic(const int &n, const double &l, const double &r, Function & func):
        ppForm_cubic(n, l, r, func, "natural") {}
};

#endif