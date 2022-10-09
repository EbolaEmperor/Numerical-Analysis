#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <iostream>
#include <vector>
#include <limits>
#include <string>
#include <sstream>

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

    NewtonInterpolation(const NewtonInterpolation & p):
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

    static const int OUTPUT_NORMAL = 0;
    static const int OUTPUT_LATEX = 1;
    static const int OUTPUT_TIKZ = 2;
    static int outputMode;

    static void setOutput(const int & style){
        if(style == OUTPUT_NORMAL){
            outputMode = OUTPUT_NORMAL;
            std::cerr << "[Interpolation output mode : Normal]" << std::endl;
        }
        else if(style == OUTPUT_LATEX){
            outputMode = OUTPUT_LATEX;
            std::cerr << "[Interpolation output mode : Latex]" << std::endl;
        }
        else if(style == OUTPUT_TIKZ){
            outputMode = OUTPUT_TIKZ;
            std::cerr << "[Interpolation output mode : Tikz]" << std::endl;
        }
        else
            std::cerr << "[Error] Incorrect interpolation output mode setting !!!" << std::endl;
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

#endif