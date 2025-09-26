#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <iostream>
#include <vector>
#include <string>
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

template<class T> int T_Polynomial<T>::outputMode = 0;

template<class T> class T_NewtonPolynomial{
private:
    static int outputMode;

protected:
    std::vector<T> x;
    std::vector<T> coef;
    int n;

public:
    T_NewtonPolynomial(){n = 0;}
    T_NewtonPolynomial(const T_NewtonPolynomial & p):
        x(p.x), coef(p.coef), n(p.n) {}
    ~T_NewtonPolynomial(){
        x.clear();
        coef.clear();
    }

    // 设置了一些输出格式，默认为直接阅读的格式Normal，也支持适合Latex公式显示的格式，以及Tikz可识别的格式
    static const int OUTPUT_NORMAL = 0;
    static const int OUTPUT_LATEX = 1;
    static const int OUTPUT_TIKZ = 2;

    static void setOutput(const int & style){
        if(style == OUTPUT_NORMAL){
            outputMode = OUTPUT_NORMAL;
            std::cerr << "[Newton T_Polynomial output mode : Normal]" << std::endl;
        }
        else if(style == OUTPUT_LATEX){
            outputMode = OUTPUT_LATEX;
            std::cerr << "[Newton T_Polynomial output mode : Latex]" << std::endl;
        }
        else if(style == OUTPUT_TIKZ){
            outputMode = OUTPUT_TIKZ;
            std::cerr << "[Newton T_Polynomial output mode : Tikz]" << std::endl;
        }
        else{
            std::cerr << "[Error] Incorrect Newton T_Polynomial output mode setting !!!" << std::endl;
            exit(-1);
        }
    }

    friend std::ostream & operator << (std::ostream & out, const T_NewtonPolynomial & ni){
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
                curs << "\\pi_{" << i+1 << "}(x)";
            }
        }
        return out;
    }

    T operator () (const T &vx) const{
        T prod = 1, ans = 0;
        for(int i = 0; i <= n; i++){
            ans += coef[i] * prod;
            prod *= vx - x[i];
        }
        return ans;
    }

    std::vector<T> operator() (const std::vector<T> &vx) const{
        int m = vx.size();
        std::vector<T> ans(m);
        for(int i = 0; i < m; i++){
            ans[i] = operator()(vx[i]);
        }
        return ans;
    }

    T_Polynomial<T> standardize(){
        T_Polynomial<T> res(coef[0]), prod(1);
        for(int i = 1; i <= n; i++){
            T_Polynomial<T> tmp(-x[i-1], 1);
            prod = prod * tmp;
            res = res + prod * coef[i];
        }
        return res;
    }
};

template<class T> int T_NewtonPolynomial<T>::outputMode = 0;

template<class T> class T_HermiteInterpolation : public T_NewtonPolynomial<T>{
private:
    std::vector<T> diffTable;
    using T_NewtonPolynomial<T>::n;
    using T_NewtonPolynomial<T>::x;
    using T_NewtonPolynomial<T>::coef;
    // 非模板类不需要上述using，但是在模板类的继承机制下，不加上面三句using会编译失败

public:
    T_HermiteInterpolation(const T_HermiteInterpolation & rhs):
        T_NewtonPolynomial<T>(rhs), diffTable(rhs.diffTable) {}
    ~T_HermiteInterpolation(){
        diffTable.clear();
    }

    // 添加插值点，若同一个插值点newx被连续添加多次，则第k个对应的newf值，会被视为newx处的(k-1)阶导数值
    void addPoint(const T &newx, const T &newf){
        n++;
        static T fac = 1;
        static int k = 0;
        if(!x.empty() && newx == x.back()){
            fac *= (++k);
        } else {
            fac = 1;
            k = 0;
        }
        x.push_back(newx);
        T newv = newf / fac;
        for(int j = 1 + k; j <= n; j++){
            if(x[n] == x[n-j]){
                std::cerr << "[Error] The interpolating point " << x[n] << " appears in a wrong position!" << std::endl;
                exit(-1);
            }
            T tmp = (newv - diffTable[j-1]) / (x[n] - x[n-j]);
            diffTable[j-1] = newv;
            newv = tmp;
        }
        diffTable.push_back(newv);
        coef.push_back(newv);
    }

    T_HermiteInterpolation(const std::vector<T> & _x, const std::vector<T> & _f){ 
        n = -1;
        if(_x.size() != _f.size()){
            std::cerr << "[Error] The size of interpolating points and interpolating values must coincide !!!" << std::endl;
            exit(-1);
        }
        int k = 0;
        for(int i = 0; i < _x.size(); i++){
            if(i && _x[i] == _x[i-1]) k++;
            else k = 0;
            addPoint(_x[i], _f[i]);
        }
    }
};

typedef T_Polynomial<double> Polynomial;
typedef T_NewtonPolynomial<double> NewtonPolynomial;
typedef T_HermiteInterpolation<double> HermiteInterpolation;
typedef HermiteInterpolation NewtonInterpolation;

#endif