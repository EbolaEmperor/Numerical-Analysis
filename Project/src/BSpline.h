#ifndef _BSPLINE_H_
#define _BSPLINE_H_

#include <iostream>
#include <vector>
#include <cstring>
#include "matrix.h"
#include "Function.h"
using std::vector;

template<int order>
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
    double d3B(const int &i, const int &k, const double &x) const{
        return k*d2B(i, k-1, x) / (knots[i+k-1] - knots[i-1]) - k*d2B(i+1, k-1, x) / (knots[i+k] - knots[i]);
    }

    double B(const int &i, const double &x) const{
        return B(i, order, x);
    };
    double dB(const int &i, const double &x) const{
        return dB(i, order, x);
    };
    double d2B(const int &i, const double &x) const{
        return d2B(i, order, x);
    };
    double d3B(const int &i, const double &x) const{
        return d2B(i, order, x);
    };

public:
    BSpline_base(){}

    BSpline_base(const BSpline_base & rhs):
        coef(rhs.coef), knots(rhs.knots) {}

    double operator () (const double &x) const{
        double ans = 0.0;
        for(int i = 0; i < coef.size(); i++)
            ans += coef[i] * B(i+1, x);
        return ans;
    }
};

class BSpline_linear : public BSpline_base<1>{
public:
    BSpline_linear(const vector<double> & t, const vector<double> & f){
        knots.clear();
        knots.push_back(t.front()-1.0);
        for(auto & x : t) knots.push_back(x);
        knots.push_back(t.back()+1.0);
        coef = f;
    }

    BSpline_linear(const vector<double> &x, Function & func){
        vector<double> f;
        for(auto & xv : x) f.push_back(func(xv));
        (*this) = BSpline_linear(x, f);
    }

    BSpline_linear(const int &n, const double &l, const double &r, Function &func){
        vector<double> t(n);
        for(int i = 0; i < n-1; i++)
            t[i] = l + (r-l)*i/(n-1);
        t.back() = r;
        (*this) = BSpline_linear(t, func);
    }
};

class BSpline_quadratic : public BSpline_base<2>{
private:
    using BSpline_base<2>::B;
    using BSpline_base<2>::dB;
    using BSpline_base<2>::d2B;

public:
    BSpline_quadratic(const vector<double> & t, const vector<double> & f){
        for(int i = 0; i < t.size()-1; i++)
            if(t[i] >= t[i+1]) error("BSpline_quadratic :: the knots must be strictly increasing!");
        
        knots.clear();
        knots.push_back(t.front()-2.0);
        knots.push_back(t.front()-1.0);
        for(auto & x : t) knots.push_back(x);
        knots.push_back(t.back()+1.0);
        knots.push_back(t.back()+2.0);

        if(f.size() != t.size()+1)
            error("BSpline_quadradic :: Wrong size of vector<double> f!");
        Matrix A(t.size()+1);
        ColVector b(t.size()+1);
        A[0][0] = B(1, t.front());
        A[0][1] = B(2, t.front());
        b[0] = f.front();
        for(int i = 0; i < t.size()-1; i++){
            double m = (t[i] + t[i+1]) / 2;
            A[i+1][i] = B(i+1, m);
            A[i+1][i+1] = B(i+2, m);
            A[i+1][i+2] = B(i+3, m);
            b[i+1] = f[i+1];
        }
        A[t.size()][t.size()-1] = B(t.size(), t.back());
        A[t.size()][t.size()] = B(t.size()+1, t.back());
        b[t.size()] = f.back();

        ColVector res = solve(A, b);
        coef.clear();
        for(int i = 0; i <= t.size(); i++)
            coef.push_back(res[i]);
    }

    BSpline_quadratic(const vector<double> &x, Function & func){
        vector<double> f;
        f.push_back(func(x.front()));
        for(int i = 0; i < x.size()-1; i++)
            f.push_back( func( (x[i]+x[i+1])/2 ) );
        f.push_back(func(x.back()));
        (*this) = BSpline_quadratic(x, f);
    }

    BSpline_quadratic(const int &n, const double &l, const double &r, Function &func){
        vector<double> t(n);
        for(int i = 0; i < n-1; i++)
            t[i] = l + (r-l)*i/(n-1);
        t.back() = r;
        (*this) = BSpline_quadratic(t, func);
    }
};

class BSpline_cubic : public BSpline_base<3>{
private:
    using BSpline_base<3>::B;
    using BSpline_base<3>::dB;
    using BSpline_base<3>::d2B;

public:
    BSpline_cubic(){}

    BSpline_cubic(const vector<double> & t, const vector<double> & f, const std::string &bondary){
        for(int i = 0; i < t.size()-1; i++)
            if(t[i] >= t[i+1]) error("BSpline_cubic :: the knots must be strictly increasing!");
        
        knots.clear();
        knots.push_back(t.front()-3.0);
        knots.push_back(t.front()-2.0);
        knots.push_back(t.front()-1.0);
        for(auto & x : t) knots.push_back(x);
        knots.push_back(t.back()+1.0);
        knots.push_back(t.back()+2.0);
        knots.push_back(t.back()+3.0);

        if(f.size() < t.size())
            error("BSpline_cubic :: Too small size of vector<double> f!");
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
        } else if(bondary == "not-a-knot"){
            if(t.size() < 4)
                error("BSpline_cubic :: not-a-knot condition needs at least 4 knots!");
            for(int j = 0; j < 3; j++){
                A[t.size()][j] = d3B(1+j, t[1]);
                A[t.size()+1][t.size()-1+j] = d2B(t.size()+j, t[t.size()-2]);
            }
            b[t.size()] = 0;
            b[t.size()+1] = 0;
        } else if(bondary == "complete"){
            if(f.size() < t.size()+2)
                error("BSpline_cubic :: Cannot read derivatives at ends in vector<double> f!");
            for(int j = 0; j < 3; j++){
                A[t.size()][j] = dB(1+j, t.front());
                A[t.size()+1][t.size()-1+j] = dB(t.size()+j, t.back());
            }
            b[t.size()] = f[t.size()];
            b[t.size()+1] = f[t.size()+1];
        } else if(bondary == "second-derivatives-at-end"){
            if(f.size() < t.size()+2)
                error("BSpline_cubic :: Cannot read second-derivatives at ends in vector<double> f!");
            for(int j = 0; j < 3; j++){
                A[t.size()][j] = d2B(1+j, t.front());
                A[t.size()+1][t.size()-1+j] = d2B(t.size()+j, t.back());
            }
            b[t.size()] = f[t.size()];
            b[t.size()+1] = f[t.size()+1];
        } else if(bondary == "periodic"){
            if(fabs(f.front() - f.back()) > 1e-14)
                std::cerr << "[Warning] BSpline_cubic :: The gap between left value and right value are larger than 1e-14 in boundary 'periodic'!" << std::endl;
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
        for(int i = 0; i < n-1; i++)
            t[i] = l + (r-l)*i/(n-1);
        t.back() = r;
        (*this) = BSpline_cubic(t, func, bondary);
    }

    BSpline_cubic(const int &n, const double &l, const double &r, Function & func):
        BSpline_cubic(n, l, r, func, "natural") {}
};

struct Point{
    double x, y;
};
double dist(const Point &a, const Point &b){
    return sqrt( (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) );
}

class Curve{
private:
    int n;
    BSpline_cubic fx, fy;

public:
    Curve() : n(0) {}
    Curve(const vector<double> &x, const vector<double> &y, const std::string & bondary){
        if(x.size() != y.size()) error("Curve :: vector x and y must be the same length!");
        n = x.size();
        vector<double> t;
        for(int i = 0; i < x.size(); i++)
            t.push_back(i);
        fx = BSpline_cubic(t, x, bondary);
        fy = BSpline_cubic(t, y, bondary);
    }

    Point operator () (const double &t){
        if(t < 0 || t > 1) error("Curve :: Out of range [0,1]!");
        if(n == 0) error("Curve :: Visited operator () at an empty curve!");
        return (Point){fx(t*(n-1)), fy(t*(n-1))};
    }

    friend std::ostream & operator << (std::ostream & out, const Curve &rhs){
        static int cnt = 0;
        cnt++;
        out << "x" << cnt << " = [ ";
        for(double i = 0; i <= rhs.n - 1 + 1e-3; i += 0.01)
            out << rhs.fx(i) << " ";
        out << "];" << std::endl;
        out << "y" << cnt << " = [ ";
        for(double i = 0; i <= rhs.n - 1 + 1e-3; i += 0.01)
            out << rhs.fy(i) << " ";
        out << "];" << std::endl;
        return out;
    }
};

#endif