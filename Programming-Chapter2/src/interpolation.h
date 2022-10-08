#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <iostream>
#include <vector>
#include <limits>

const double _epsL = 10 * std::numeric_limits<double>::epsilon();

class Function{
public:
    virtual double operator () (const double &x) const = 0;
    virtual double diff (const double &x) const{
        return ((*this)(x+_epsL)-(*this)(x-_epsL)) / (2*_epsL);
    }
};

class NewtonInterpolation{
private:
    Function & f;
    std::vector<double> x;
    std::vector<double> diffTable;
    std::vector<double> coef;
    int n;

public:
    void addPoint(const double &newx){
        n++;
        x.push_back(newx);
        double newv = f(newx);
        for(int j = 1; j <= n; j++){
            double tmp = (newv - diffTable[j-1]) / (x[n] - x[n-j]);
            diffTable[j-1] = newv;
            newv = tmp;
        }
        diffTable.push_back(newv);
        coef.push_back(newv);
    }

    NewtonInterpolation(NewtonInterpolation & p):
        f(p.f), x(p.x), diffTable(p.diffTable), coef(p.coef), n(p.n) {}
    NewtonInterpolation(Function & f, std::vector<double> & x): f(f) { 
        n = -1;
        for(double v : x) addPoint(v);
    }
    ~NewtonInterpolation(){
        x.clear();
        diffTable.clear();
        coef.clear();
    }

    double operator () (const double &vx) const{
        double prod = 1, ans = 0;
        for(int i = 0; i <= n; i++){
            ans += coef[i] * prod;
            prod *= vx - x[i];
        }
        return ans;
    }
};

#endif