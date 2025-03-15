#ifndef _TRIANGULAR_INTERPOLATION_H_
#define _TRIANGULAR_INTERPOLATION_H_

#include <vector>
#include <complex>
#include <string>

class Vector{
private:
    std::vector<double> values;
public:
    int size() const;
    double& operator [] (const int &i);
    const double& operator [] (const int &i) const;
    void push_back(const double &x);
    void load(const std::string& fname);
};

class TriangularInterpolator{
private:
    std::vector< std::complex<double> > coef;
    std::vector< std::complex<double> > fft(const std::vector< std::complex<double> > &, const int &v) const;
public:
    double at(const double &t) const;
    void fit(const Vector &x);
    void filter(const int &cutoff);
    Vector recover() const;
    void load(const std::string& fname);
    void save(const std::string& fname) const;
    void discreteOutput(const std::string& fname, const int &p) const;
};

#endif