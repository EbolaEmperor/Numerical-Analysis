#include <iostream>
#include <fstream>
#include <cmath>
#include "ArnoldiInterp.h"
#include "interpolation.h"
#include "simpson.h"
#include "plot.h"
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
    std::vector<double> x(n), f(n);
    for(int i = 1; i <= n; i++){
        x[i-1] = cos((2.0 * i - 1.0) / (2.0 * n) * M_PI);
        f[i-1] = func(x[i-1]);
    }
    auto poly = ArnoldiInterp(x, f);
    auto npoly = NewtonInterpolation(x, f);

    const int N = 100;
    std::vector<double> grid(2 * N + 1);
    for (int i = -N; i <= N; i++)
        grid[i + N] = 1.0 * i / N;
    auto val = poly(grid);

    plotcpp::Plot plot(true);
    // plot.SetTerminalRaw("pngcairo size 1600,1200");
    // plot.SetOutput("arnoldi.png");
    plot.SetTerminal("qt");
    plot.SetXLabel("x");
    plot.SetYLabel("y");
    plot.SetXRange(-1, 1);
    plot.SetYRange(-0.3, 1);

    auto y_arnoldi = poly(grid);
    auto y_newton = npoly(grid);

    plot.Draw2D(plotcpp::Lines(grid.begin(), grid.end(), y_newton.begin(), "Newton", plotcpp::LineType::Solid, 3),
                plotcpp::Lines(grid.begin(), grid.end(), y_arnoldi.begin(), "Arnoldi", plotcpp::LineType::DashDot, 3));
    plot.Flush();

    auto L2errfun = [&](double x) {
        double tmp = func(x) - poly(x);
        return tmp * tmp;
    };
    double L2err = adpt_simpson(L2errfun, -1.0, 1.0, 2.3e-16);
    cout << "L2 error (Arnoldi) = " << sqrt(L2err) << endl;

    return 0;
}