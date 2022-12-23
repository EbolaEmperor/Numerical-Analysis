#ifndef _BSPLINE_H_
#define _BSPLINE_H_

#include <iostream>
#include <vector>
#include <cstring>
#include "matrix.h"
#include "Function.h"
using std::vector;

class BSpline_base{
protected:
    vector<double> coef;
    vector<double> knots;
    double B(const int &i, const int &k, const double &x) const{
        return k==0 ? (knots[i-1]<x && x<=knots[i]) : 
            (x-knots[i-1])/(knots[i+k-1]-knots[i-1])*B(i,k-1,x) + (knots[i+k]-x)/(knots[i+k]-knots[i])*B(i+1,k-1,x);
    }
    double dB(const int &i, const int &k, const double &x) const{
        return k*B(i, k-1, x) / (knots[i+k-1] - knots[i-1]) - k*B(i+1, k-1, x) / (knots[i+k] - knots[i]);
    }
    double d2B(const int &i, const int &k, const double &x) const{
        return k*dB(i, k-1, x) / (knots[i+k-1] - knots[i-1]) - k*dB(i+1, k-1, x) / (knots[i+k] - knots[i]);
    }

    virtual double B(const int &i, const double &x) const = 0;
    virtual double dB(const int &i, const double &x) const = 0;
    virtual double d2B(const int &i, const double &x) const = 0;

public:
    BSpline_base(){}
    BSpline_base(const BSpline_base & rhs):
        coef(rhs.coef), knots(rhs.knots) {}

    double operator () (const double &x){
        double ans = 0.0;
        for(int i = 0; i < coef.size(); i++)
            ans += coef[i] * B(i+1, x);
        return ans;
    }
};

class BSpline_linear : public BSpline_base{
private:
    using BSpline_base::B;
    double B(const int &i, const double &x) const{
        return B(i, 1, x);
    }

public:
    BSpline_linear(const vector<double> & t, const vector<double> & f){
        knots.clear();
        knots.push_back(t.front()-1.0);
        for(auto & x : t) knots.push_back(x);
        knots.push_back(t.back()+1.0);
        coef = f;
    }
};

class BSpline_cubic : public BSpline_base{
private:
    using BSpline_base::B;
    using BSpline_base::dB;
    using BSpline_base::d2B;
    
    double B(const int &i, const double &x) const{
        return B(i, 3, x);
    }
    double dB(const int &i, const double &x) const{
        return dB(i, 3, x);
    }
    double d2B(const int &i, const double &x) const{
        return d2B(i, 3, x);
    }

public:
    BSpline_cubic(const vector<double> & t, const vector<double> & f, const std::string &bondary){
        knots.clear();
        knots.push_back(t.front()-3.0);
        knots.push_back(t.front()-2.0);
        knots.push_back(t.front()-1.0);
        for(auto & x : t) knots.push_back(x);
        knots.push_back(t.back()+1.0);
        knots.push_back(t.back()+2.0);
        knots.push_back(t.back()+3.0);

        Matrix A(t.size()+2);
        ColVector b(t.size()+2);
        for(int i = 0; i < t.size(); i++){
            for(int j = 0; j < 3; j++)
                A[i][i+j] = B(i+1+j, t[i]);
            b[i] = f[i];
        }

        if(bondary == "natural"){
            for(int j = 0; j < 3; j++){
                A[t.size()][j] = d2B(1+j, t.front());
                A[t.size()+1][t.size()-1+j] = d2B(t.size()+j, t.back());
            }
            b[t.size()] = 0;
            b[t.size()+1] = 0;
        } else if(bondary == "complete"){
            for(int j = 0; j < 3; j++){
                A[t.size()][j] = dB(1+j, t.front());
                A[t.size()+1][t.size()-1+j] = dB(t.size()+j, t.back());
            }
            b[t.size()] = f[t.size()];
            b[t.size()+1] = f[t.size()+1];
        } else if(bondary == "second-derivatives-at-end"){
            for(int j = 0; j < 3; j++){
                A[t.size()][j] = d2B(1+j, t.front());
                A[t.size()+1][t.size()-1+j] = d2B(t.size()+j, t.back());
            }
            b[t.size()] = f[t.size()];
            b[t.size()+1] = f[t.size()+1];
        } else if(bondary == "periodic"){
            for(int j = 0; j < 3; j++){
                A[t.size()][j] = dB(1+j, t.front());
                A[t.size()][t.size()-1+j] = -dB(t.size()+j, t.back());
                A[t.size()+1][j] = d2B(1+j, t.front());
                A[t.size()+1][t.size()-1+j] = -d2B(t.size()+j, t.back());
            }
            b[t.size()] = b[t.size()+1] = 0;
        }

        ColVector res = solve(A, b);
        coef.clear();
        for(int i = 0; i < t.size()+2; i++)
            coef.push_back(res[i]);
    }

    BSpline_cubic(const vector<double> & t, const vector<double> & f):
        BSpline_cubic(t, f, "natural") {}

    BSpline_cubic(const vector<double> & t, Function & func, const std::string &bondary){
        vector<double> f;
        for(auto & x : t) f.push_back(func(x));
        if(bondary == "complete"){
            f.push_back(func.diff(t.front()));
            f.push_back(func.diff(t.back()));
        } else if(bondary == "second-derivatives-at-end"){
            f.push_back(func.diff2(t.front()));
            f.push_back(func.diff2(t.back()));
        }
        (*this) = BSpline_cubic(t, f, bondary);
    }

    BSpline_cubic(const vector<double> & t, Function & func):
        BSpline_cubic(t, func, "natural") {}

    BSpline_cubic(const int &n, const double &l, const double &r, Function & func, const std::string &bondary){
        vector<double> t(n);
        for(int i = 0; i < n; i++)
            t[i] = l + (r-l)*i/(n-1);
        (*this) = BSpline_cubic(t, func, bondary);
    }

    BSpline_cubic(const int &n, const double &l, const double &r, Function & func):
        BSpline_cubic(n, l, r, func, "natural") {}
};

#endif