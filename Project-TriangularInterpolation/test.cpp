#include "TriangularInterpolation.h"
#include <iostream>
#include <fstream>
using namespace std;

int main(){
    Vector x;
    x.load("data/input.in");
    TriangularInterpolator func;
    func.fit(x);
    func.save("data/sample.model");
    func.discreteOutput("data/sample.out", 256);

    Vector rx = func.recover();
    cout << "Recover test: " << endl;
    for(int i = 0; i < rx.size(); i++)
        cout << rx[i] << " ";
    cout << endl << endl;

    TriangularInterpolator refunc;
    refunc.load("data/sample.model");
    rx = refunc.recover();
    cout << "Load and recover test: " << endl;
    for(int i = 0; i < rx.size(); i++)
        cout << rx[i] << " ";
    cout << endl;
    cout << endl << "Run the following command in gnuplot to see the fitted curve:" << endl << endl;
    cout << "   plot \"data/sample.out\" w l title \"fitted curve\",\\" << endl
         << "        \"data/input.gnuplot_in\" w p pt 6 title \"input data\"" << endl << endl;
    return 0;
}