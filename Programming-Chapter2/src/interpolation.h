#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "matrix.h"

class NewtonInterpolation{
private:
    std::vector<double> x;
    std::vector<double> diffTable;
    std::vector<double> coef;
    int n;

public:
    // 添加一个插值点，复杂度为O(n)，允许在构造函数外调用
    void addPoint(const double &newx, const double &newf){
        n++;
        x.push_back(newx);
        double newv = newf;
        for(int j = 1; j <= n; j++){
            double tmp = (newv - diffTable[j-1]) / (x[n] - x[n-j]);
            diffTable[j-1] = newv;
            newv = tmp;
        }
        diffTable.push_back(newv);
        coef.push_back(newv);
    }

    NewtonInterpolation(const NewtonInterpolation & p):
        x(p.x), diffTable(p.diffTable), coef(p.coef), n(p.n) {}
    NewtonInterpolation(std::vector<double> & _x, std::vector<double> & _f){ 
        n = -1;
        if(_x.size() != _f.size()){
            std::cerr << "[Error] The size of interpolating points and interpolating values must coincide !!!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i < _x.size(); i++)
            addPoint(_x[i], _f[i]);
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

    // 设置了一些输出格式，默认为直接阅读的格式Normal，也支持适合Latex公式显示的格式，以及Tikz可识别的格式
    static const int OUTPUT_NORMAL = 0;
    static const int OUTPUT_LATEX = 1;
    static const int OUTPUT_TIKZ = 2;
    static int outputMode;

    static void setOutput(const int & style){
        if(style == OUTPUT_NORMAL){
            outputMode = OUTPUT_NORMAL;
            std::cerr << "[Newton interpolation output mode : Normal]" << std::endl;
        }
        else if(style == OUTPUT_LATEX){
            outputMode = OUTPUT_LATEX;
            std::cerr << "[Newton interpolation output mode : Latex]" << std::endl;
        }
        else if(style == OUTPUT_TIKZ){
            outputMode = OUTPUT_TIKZ;
            std::cerr << "[Newton Interpolation output mode : Tikz]" << std::endl;
        }
        else{
            std::cerr << "[Error] Incorrect Newton interpolation output mode setting !!!" << std::endl;
            exit(-1);
        }
    }

    friend std::ostream & operator << (std::ostream & out, NewtonInterpolation & ni){
        std::stringstream curs;
        for(int i = 0; i <= ni.n; i++){
            out << ni.coef[i] << curs.str();
            if(i < ni.n && ni.coef[i+1] >= 0) out << "+";
            if(outputMode == OUTPUT_TIKZ){
                if(ni.x[i] == 0)
                    curs << "*\\x";
                else if(ni.x[i] > 0)
                    curs << "*(\\x-" << ni.x[i] << ")";
                else
                    curs << "*(\\x+" << -ni.x[i] << ")";
            }
            else if(outputMode == OUTPUT_NORMAL){
                if(ni.x[i] == 0)
                    curs << "*x";
                else if(ni.x[i] > 0)
                    curs << "*(x-" << ni.x[i] << ")";
                else
                    curs << "*(x+" << -ni.x[i] << ")";
            }
            else if(outputMode == OUTPUT_LATEX){
                curs.str("");
                curs.clear();
                curs << "\\pi_{" << i << "}(x)";
            }
        }
        return out;
    }
};

int NewtonInterpolation::outputMode = 0;

class HermiteInterpolation{
private:
    std::vector<double> coef;
    int n;
    
    struct Condition{
        double x;
        int order;
        double value;
    };
    std::vector<Condition> cond;

    void solve(){
        if(!coef.empty()) coef.clear();
        n = cond.size() - 1;
        Matrix A(n+1, n+1);
        ColVector b(n+1);
        for(int i = 0; i <= n; i++){
            double x = cond[i].x, prod = 1.0;
            int k = cond[i].order;
            b[i] = cond[i].value;
            for(int j = k; j <= n; j++){
                double fac = 1.0;
                for(int p = 0; p < k; p++)
                    fac *= j - p;
                A[i][j] = prod * fac;
                prod *= x;
            }
        }
        ColVector res = ::solve(A, b);
        if(res.empty()){
            std::cerr << "[Error] Cannot solve Hermite interpolation, check your condition!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i <= n; i++)
            coef.push_back(res[i]);
    }

public:
    HermiteInterpolation(): n(0) {}
    HermiteInterpolation(const HermiteInterpolation & rhs):
        cond(rhs.cond), coef(rhs.coef), n(rhs.n) {}
    
    void addCondition(const double &x, const int &order, const double &value){
        cond.push_back((Condition){x,order,value});
        if(!coef.empty()) coef.clear();
    }
    
    double operator () (const double &x) {
        if(coef.empty()) solve();
        double prod = 1, ans = 0;
        for(int i = 0; i <= n; i++){
            ans += coef[i] * prod;
            prod *= x;
        }
        return ans;
    }

    // 设置了一些输出格式，支持适合Latex公式显示的格式，以及Tikz可识别的格式
    static const int OUTPUT_LATEX = 0;
    static const int OUTPUT_TIKZ = 1;
    static int outputMode;

    static void setOutput(const int & style){
        if(style == OUTPUT_LATEX){
            outputMode = OUTPUT_LATEX;
            std::cerr << "[Hermite interpolation output mode : Latex]" << std::endl;
        }
        else if(style == OUTPUT_TIKZ){
            outputMode = OUTPUT_TIKZ;
            std::cerr << "[Hermite interpolation output mode : Tikz]" << std::endl;
        }
        else{
            std::cerr << "[Error] Incorrect Hermite interpolation output mode setting !!!" << std::endl;
            exit(-1);
        }
    }

    friend std::ostream & operator << (std::ostream & out, HermiteInterpolation & ni){
        if(ni.coef.empty()) ni.solve();
        std::stringstream curs;
        for(int i = 0; i <= ni.n; i++){
            out << ni.coef[i] << curs.str();
            if(i < ni.n && ni.coef[i+1] >= 0) out << "+";
            if(outputMode == OUTPUT_TIKZ)
                curs << "*\\x";
            else if(outputMode == OUTPUT_LATEX){
                curs.str("");
                curs.clear();
                curs << "x^{" << i + 1 << "}";
            }
        }
        return out;
    }
};

int HermiteInterpolation::outputMode = 0;

#endif