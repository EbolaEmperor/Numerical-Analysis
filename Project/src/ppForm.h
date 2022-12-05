#ifndef _PPFORM_H_
#define _PPFORM_H_

#include <vector>
#include <iostream>
#include <cstring>
#include <cmath>
#include <limits>
#include "Polynimial.h"
#include "matrix.h"

const double _epsL = 10 * std::numeric_limits<double>::epsilon();

// 函数虚类，实体函数需要继承此类，并定义()运算，子类可以定义diff为导函数，若不定义，则diff默认采用差商代替导数的方法
class Function{
public:
    virtual double operator () (const double &x) const = 0;

    virtual double diff (const double &x) const{
        return ((*this)(x+_epsL)-(*this)(x-_epsL)) / (2*_epsL);
    }

    virtual double diff2 (const double &x) const{
        return ((*this)(x+2*_epsL)+(*this)(x-2*_epsL)-2*(*this)(x)) / (4*_epsL*_epsL);
    }
};

class ppForm_base{
protected:
    int n;
    std::vector<double> knots;
    std::vector<Polynomial> poly;

public:
    ppForm_base() {}

    ppForm_base(const int &n): n(n){
        knots.resize(n+1);
        poly.resize(n);
    }

    ppForm_base(const ppForm_base &rhs):
        n(rhs.n), knots(rhs.knots), poly(rhs.poly) {}
    
    double operator () (const double &x) const{
        if(n == 0){
            std::cerr << "[Error] ppForm operator () :: Visited an uninitialized object." << std::endl;
            exit(-1);
        }
        if(x < knots.front() || x > knots.back()){
            std::cerr << "[Error] ppForm operator () :: point x are not in the interpolation interval." << std::endl;
            exit(-1);
        }
        for(int i = 0; i < n-1; i++){
            if(x <= knots[i+1]) return poly[i](x);
        }
        return poly[n-1](x);
    }
};

class ppForm_linear : public ppForm_base{
public:
    ppForm_linear(const std::vector<double> &x, const std::vector<double> &f){
        if(x.size() != f.size() || x.size() < 2){
            std::cerr << "[Error] ppForm_linear initializing :: The knots and values are incorrect." << std::endl;
            exit(-1);
        }
        n = x.size() - 1;
        for(int i = 0; i < n; i++){
            if(x[i] >= x[i+1]){
                std::cerr << "[Error] ppForm_linear initializing :: The knots are not increasing." << std::endl;
                exit(-1);
            }
        }
        knots = x;
        for(int i = 0; i < n; i++){
            poly.push_back( Polynomial( (x[i+1]*f[i]-x[i]*f[i+1])/(x[i+1]-x[i]), (f[i+1]-f[i])/(x[i+1]-x[i]) ) );
        }
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
            for(int j = 2; j < 4; j++){
                A[4*n-2][j] = j*(j-1) * pow(x[0], j-2);
                A[4*n-1][4*(n-1)+j] = j*(j-1) * pow(x[n], j-2);
            }
            b[4*n-2] = f[n+1];
            b[4*n-1] = f[n+2];
        } else if(bondary == "periodic"){
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
        //std::cerr << A << std::endl << b.T() << std::endl;
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
};

#endif