#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

class Polynomial{
private:
    static int outputMode;

protected:
    std::vector<double> coef;
    int n;

public:
    Polynomial(){n = 0;}
    Polynomial(const double &x){
        n = 0;
        coef.push_back(x);
    }
    Polynomial(const double &x0, const double &x1){
        n = 1;
        coef.push_back(x0);
        coef.push_back(x1);
    }
    Polynomial(const Polynomial & p):
        coef(p.coef), n(p.n) {}

    // 设置了一些输出格式，默认为适合Latex公式显示的格式，也支持Tikz可识别的格式
    static const int OUTPUT_LATEX = 0;
    static const int OUTPUT_TIKZ = 1;

    static void setOutput(const int & style){
        if(style == OUTPUT_LATEX){
            outputMode = OUTPUT_LATEX;
            std::cerr << "[Polynomial output mode : Latex]" << std::endl;
        }
        else if(style == OUTPUT_TIKZ){
            outputMode = OUTPUT_TIKZ;
            std::cerr << "[Polynomial output mode : Tikz]" << std::endl;
        }
        else{
            std::cerr << "[Error] Incorrect Polynomial output mode setting !!!" << std::endl;
            exit(-1);
        }
    }

    friend std::ostream & operator << (std::ostream & out, Polynomial & ni){
        std::stringstream curs;
        for(int i = 0; i <= ni.n; i++){
            out << ni.coef[i] << curs.str();
            if(i < ni.n && ni.coef[i+1] >= 0) out << "+";
            if(outputMode == OUTPUT_TIKZ)
                curs << "*\\x";
            else if(outputMode == OUTPUT_LATEX){
                curs.str("");
                curs.clear();
                curs << "*x^{" << i + 1 << "}";
            }
        }
        return out;
    }

    double operator () (const double &vx) const{
        double prod = 1, ans = 0;
        for(int i = 0; i <= n; i++){
            ans += coef[i] * prod;
            prod *= vx;
        }
        return ans;
    }

    Polynomial operator + (const Polynomial &rhs){
        Polynomial res(*this);
        res.n = std::max(res.n, rhs.n);
        res.coef.resize(res.n + 1);
        for(int i = 0; i <= rhs.n; i++)
            res.coef[i] += rhs.coef[i];
        return res;
    }

    Polynomial operator * (const double &rhs){
        Polynomial res;
        res.n = n;
        res.coef.resize(res.n + 1);
        for(int i = 0; i <= n; i++)
            res.coef[i] = coef[i] * rhs;
        return res;
    }

    Polynomial operator * (const Polynomial &rhs){
        Polynomial res;
        res.n = n + rhs.n;
        res.coef.resize(res.n + 1);
        for(int i = 0; i <= n; i++)
            for(int j = 0; j <= rhs.n; j++)
                res.coef[i+j] += coef[i] * rhs.coef[j];
        return res;
    }

    Polynomial diff(){
        Polynomial res;
        res.n = n ? n - 1 : 0;
        for(int i = 1; i <= n; i++)
            res.coef.push_back(coef[i] * i);
        return res;
    }
};

int Polynomial::outputMode = 0;

class NewtonPolynomial{
private:
    static int outputMode;

protected:
    std::vector<double> x;
    std::vector<double> diffTable;
    std::vector<double> coef;
    int n;

public:
    NewtonPolynomial(){n = 0;}
    NewtonPolynomial(const NewtonPolynomial & p):
        x(p.x), diffTable(p.diffTable), coef(p.coef), n(p.n) {}
    ~NewtonPolynomial(){
        x.clear();
        diffTable.clear();
        coef.clear();
    }

    // 设置了一些输出格式，默认为直接阅读的格式Normal，也支持适合Latex公式显示的格式，以及Tikz可识别的格式
    static const int OUTPUT_NORMAL = 0;
    static const int OUTPUT_LATEX = 1;
    static const int OUTPUT_TIKZ = 2;

    static void setOutput(const int & style){
        if(style == OUTPUT_NORMAL){
            outputMode = OUTPUT_NORMAL;
            std::cerr << "[Newton polynomial output mode : Normal]" << std::endl;
        }
        else if(style == OUTPUT_LATEX){
            outputMode = OUTPUT_LATEX;
            std::cerr << "[Newton polynomial output mode : Latex]" << std::endl;
        }
        else if(style == OUTPUT_TIKZ){
            outputMode = OUTPUT_TIKZ;
            std::cerr << "[Newton polynomial output mode : Tikz]" << std::endl;
        }
        else{
            std::cerr << "[Error] Incorrect Newton polynomial output mode setting !!!" << std::endl;
            exit(-1);
        }
    }

    friend std::ostream & operator << (std::ostream & out, NewtonPolynomial & ni){
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

    double operator () (const double &vx) const{
        double prod = 1, ans = 0;
        for(int i = 0; i <= n; i++){
            ans += coef[i] * prod;
            prod *= vx - x[i];
        }
        return ans;
    }

    Polynomial standardize(){
        Polynomial res(coef[0]), prod(1.0);
        for(int i = 1; i <= n; i++){
            Polynomial tmp(-x[i-1], 1.0);
            prod = prod * tmp;
            res = res + prod * coef[i];
        }
        return res;
    }
};

int NewtonPolynomial::outputMode = 0;

class NewtonInterpolation : public NewtonPolynomial{
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

    NewtonInterpolation(std::vector<double> & _x, std::vector<double> & _f){ 
        n = -1;
        if(_x.size() != _f.size()){
            std::cerr << "[Error] The size of interpolating points and interpolating values must coincide !!!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i < _x.size(); i++)
            addPoint(_x[i], _f[i]);
    }
};

class HermiteInterpolation : public NewtonPolynomial{
public:
    // 添加一个插值点，复杂度为O(n)，允许在构造函数外调用
    void addPoint(const double &newx, const double &newf, const double &newdf){
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
        n++;
        x.push_back(newx);
        newv = newdf;
        for(int j = 2; j <= n; j++){
            double tmp = (newv - diffTable[j-1]) / (x[n] - x[n-j]);
            diffTable[j-1] = newv;
            newv = tmp;
        }
        diffTable.push_back(newv);
        coef.push_back(newv);
    }

    HermiteInterpolation(std::vector<double> & _x, std::vector<double> & _f, std::vector<double> & _df){ 
        n = -1;
        if(_x.size() != _f.size() || _x.size() != _df.size()){
            std::cerr << "[Error] The size of interpolating points and interpolating values must coincide !!!" << std::endl;
            exit(-1);
        }
        for(int i = 0; i < _x.size(); i++)
            addPoint(_x[i], _f[i], _df[i]);
    }
};

#endif