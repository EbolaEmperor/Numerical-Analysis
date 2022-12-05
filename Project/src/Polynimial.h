#ifndef _POLYNIMIAL_H_
#define _POLYNOMIAL_H_

#include <vector>
#include <iostream>
#include <sstream>

template<class T> class T_Polynomial{
private:
    static int outputMode;

protected:
    std::vector<T> coef;
    int n;

public:
    T_Polynomial(){n = 0;}
    T_Polynomial(const T &x){
        n = 0;
        coef.push_back(x);
    }
    T_Polynomial(const T &x0, const T &x1){
        n = 1;
        coef.push_back(x0);
        coef.push_back(x1);
    }
    T_Polynomial(const std::vector<T> & coef):
        coef(coef), n(coef.size()-1) {}
    T_Polynomial(const T_Polynomial & p):
        coef(p.coef), n(p.n) {}

    // 设置了一些输出格式，默认为适合Latex公式显示的格式，也支持Tikz可识别的格式
    static const int OUTPUT_LATEX = 0;
    static const int OUTPUT_TIKZ = 1;

    static void setOutput(const int & style){
        if(style == OUTPUT_LATEX){
            outputMode = OUTPUT_LATEX;
            std::cerr << "[T_Polynomial output mode : Latex]" << std::endl;
        }
        else if(style == OUTPUT_TIKZ){
            outputMode = OUTPUT_TIKZ;
            std::cerr << "[T_Polynomial output mode : Tikz]" << std::endl;
        }
        else{
            std::cerr << "[Error] Incorrect T_Polynomial output mode setting !!!" << std::endl;
            exit(-1);
        }
    }

    friend std::ostream & operator << (std::ostream & out, const T_Polynomial & ni){
        std::stringstream curs;
        bool first = true;
        for(int i = 0; i <= ni.n; i++){
            if(ni.coef[i] != 0){
                if(!first && ni.coef[i] >= 0) out << "+";
                out << ni.coef[i] << curs.str();
                if(first) first = false;
            }
            if(outputMode == OUTPUT_TIKZ)
                curs << "*\\x";
            else if(outputMode == OUTPUT_LATEX){
                curs.str("");
                curs.clear();
                if(i == 0) curs << "x";
                else curs << "x^{" << i + 1 << "}";
            }
        }
        return out;
    }

    T operator () (const T &vx) const{
        T prod = 1, ans = 0;
        for(int i = 0; i <= n; i++){
            ans += coef[i] * prod;
            prod *= vx;
        }
        return ans;
    }

    T_Polynomial operator + (const T_Polynomial &rhs){
        T_Polynomial res(*this);
        res.n = std::max(res.n, rhs.n);
        res.coef.resize(res.n + 1);
        for(int i = 0; i <= rhs.n; i++)
            res.coef[i] += rhs.coef[i];
        return res;
    }

    T_Polynomial operator * (const T &rhs){
        T_Polynomial res;
        res.n = n;
        res.coef.resize(res.n + 1);
        for(int i = 0; i <= n; i++)
            res.coef[i] = coef[i] * rhs;
        return res;
    }

    T_Polynomial operator * (const T_Polynomial &rhs){
        T_Polynomial res;
        res.n = n + rhs.n;
        res.coef.resize(res.n + 1);
        for(int i = 0; i <= n; i++)
            for(int j = 0; j <= rhs.n; j++)
                res.coef[i+j] += coef[i] * rhs.coef[j];
        return res;
    }

    T_Polynomial diff(){
        T_Polynomial res;
        res.n = n ? n - 1 : 0;
        for(int i = 1; i <= n; i++)
            res.coef.push_back(coef[i] * i);
        return res;
    }
};

typedef T_Polynomial<double> Polynomial;

#endif