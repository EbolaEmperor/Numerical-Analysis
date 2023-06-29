#include "TriangularInterpolation.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>

using namespace std;
typedef complex<double> Complex;

void Vector::load(const std::string& fname){
    ifstream fin(fname);
    double val;
    while(fin >> val){
        values.push_back(val);
    }
}

int Vector::size() const{
    return values.size();
}

void Vector::push_back(const double &x){
    values.push_back(x);
}

double& Vector::operator[] (const int &i){
    return values[i];
}

const double& Vector::operator[] (const int &i) const{
    return values[i];
}

vector<Complex> TriangularInterpolator::fft(const vector<Complex> &in_a, const int &v) const{
    const int n = in_a.size();
    static const double pi = acos(-1);
    int l = 0, *r = new int[n];
    memset(r, 0, sizeof(int) * n);
    for(int i = 1; i < n; i <<= 1) l++;
	for(int i = 0; i < n; i++)
        r[i] = ( r[i/2] / 2 ) | ( (i&1) << (l-1) );
    vector<Complex> a = in_a;
	for(int i = 0; i < n; i++)
        if(i < r[i]) swap(a[i], a[r[i]]);
	for(int i = 1; i < n; i <<= 1){
		Complex wn(cos(pi/i), v*sin(pi/i));
		int p = i << 1;
		for(int j = 0; j < n; j += p)
		{
			Complex w(1, 0);
			for(int k = 0; k < i; k++)
			{
				Complex x = a[j+k], y = w*a[i+j+k];
				a[j+k] = x + y;
                a[i+j+k] = x - y;
				w = w * wn;
			}
		}
	}
    return a;
}

double TriangularInterpolator::at(const double &t) const{
    const int n = coef.size();
    static const double pi = acos(-1);
    double res = 0;
    for(int k = 0; k < n; k++){
        res += coef[k].real() * cos(2*pi*k*t) - coef[k].imag() * sin(2*pi*k*t);
    }
    return res / n;
}

void TriangularInterpolator::fit(const Vector &x){
    vector<Complex> a;
    for(int i = 0; i < x.size(); i++)
        a.push_back( Complex(x[i], 0) );
    coef = fft(a, -1);
}

Vector TriangularInterpolator::recover() const{
    auto res = fft(coef, 1);
    Vector points;
    const int n = res.size();
    for(int i = 0; i < res.size(); i++)
        points.push_back(res[i].real()/n);
    return points;
}

void TriangularInterpolator::discreteOutput(const std::string& fname, const int &p) const{
    ofstream out(fname);
    vector<Complex> yp(p);
    const int n = coef.size();
    for(int i = 0; i <= n/2; i++)
        yp[i] = coef[i];
    for(int i = 0; i <= n/2-2; i++)
        yp[p-n/2+1+i] = coef[n/2+1+i];
    auto xp = fft(yp, 1);
    double delta = 1.0 / p;
    for(int i = 0; i < p; i++){
        out << i*delta << " " << xp[i].real()/n << endl;
    }
    out.close();
}

void TriangularInterpolator::save(const std::string& fname) const{
    ofstream writeFile;
    const int write_size = sizeof(double);
    writeFile.open(fname, ios::out | ios::binary);
    for (int ind = 0 ; ind < coef.size(); ind++) {
        double r = coef[ind].real();
        double i = coef[ind].imag();
        writeFile.write(reinterpret_cast<char*>(&r), write_size);
        writeFile.write(reinterpret_cast<char*>(&i), write_size);
    }
    cout << "model has been saved to " << fname << endl << endl;
}

void TriangularInterpolator::load(const std::string &fname){
    ifstream readFile;
    double temp, pre;
    bool first = true;
    coef.clear();
    readFile.open(fname, ios::in | ios::binary) ;
    while(readFile.read(reinterpret_cast<char*> (&temp), sizeof(double))) {
        if(first){
            pre = temp;
            first = false;
        } else {
            coef.push_back( Complex(pre, temp) );
            first = true;
        }
    }
}