#include "TriangularInterpolation.h"
#include <iostream>
#include <fstream>
#include <random>
using namespace std;

double func(const double &x){
    return 3 * sin(2 * M_PI * x) + cos(4 * M_PI * x) + 2 * sin(8 * M_PI * x);
}

random_device rd;
mt19937 gen(rd());
normal_distribution<double> dis(0, 0.1);

int main(){
    Vector x;
    const int N = 1024;
    for(int i = 0; i < N; i++){
        x.push_back(func(1.0 * i / N) + dis(gen));
    }

    ofstream out("data/input2.gnuplot_in");
    for(int i = 0; i < N; i++)
        out << i / double(N) << " " << x[i] << endl;
    out.close();

    TriangularInterpolator interpolator;
    interpolator.fit(x);
    interpolator.filter(4);
    interpolator.discreteOutput("data/sample.out", 2048);

    cout << endl;
    cout << endl << "Run the following command in gnuplot to see the fitted curve:" << endl << endl;
    cout << "   plot \"data/input2.gnuplot_in\" w l title \"input data\",\\" << endl
         << "        \"data/sample.out\" w l lw 2 title \"fitted curve\"" << endl << endl;
    return 0;
}