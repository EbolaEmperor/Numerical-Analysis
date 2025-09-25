#include <iostream>
#include <fstream>
#include <cmath>
#include "ArnoldiInterp.h"
#include "simpson.h"
using namespace std;
using namespace Eigen;

class Func{
public:
    double operator () (const double &x) const{
        return 1.0 / (1.0 + 25.0 * x * x);
    }
} func;

int main(){
    int n;
    cin >> n; ++n;
    VectorXd x(n), f(n);
    for(int i = 1; i <= n; i++){
        x(i-1) = cos((2.0 * i - 1.0) / (2.0 * n) * M_PI);
        f(i-1) = func(x[i-1]);
    }
    auto poly = ArnoldiInterp(x, f);

    const int N = 100;
    VectorXd grid(2 * N + 1);
    for (int i = -N; i <= N; i++)
        grid(i + N) = 1.0 * i / N;
    auto val = poly(grid);

    ofstream fout("output.txt");
    for (int i = 0; i < 2 * N + 1; i++)
        fout << grid(i) << " " << val(i) << "\n";
    fout.close();
    cout << "Plot data has been written to 'output.txt'." << endl;

    auto L2errfun = [&](double x) {
        double tmp = func(x) - poly(x);
        return tmp * tmp;
    };
    double L2err = adpt_simpson(L2errfun, -1.0, 1.0, 2.3e-16);
    cout << "L2 error = " << sqrt(L2err) << endl;

    return 0;
}