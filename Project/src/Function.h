#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include <limits>
#include <iostream>
#include <cstring>

const double _epsL = 10 * std::numeric_limits<double>::epsilon();

// 函数虚类，实体函数需要继承此类，并定义()运算，子类可以定义diff为导函数，若不定义，则diff默认采用差商代替导数的方法
class Function{
public:
    virtual double operator () (const double &x) const = 0;

    double diff (const double &x) const{
        return ((*this)(x+_epsL)-(*this)(x-_epsL)) / (2*_epsL);
    }

    double diff2 (const double &x) const{
        return ((*this)(x+2*_epsL)+(*this)(x-2*_epsL)-2*(*this)(x)) / (4*_epsL*_epsL);
    }
};

void error(const std::string & errcode){
    std::cerr << "[Error] " << errcode << std::endl;
    exit(-1);
}

#endif