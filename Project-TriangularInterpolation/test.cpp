#include "TriangularInterpolation.h"
#include <iostream>
using namespace std;

int main(){
    Vector x;
    x.load("data/input.in");
    TriangularInterpolator func;
    func.fit(x);
    func.save("data/sample.model");
    func.discreteOutput("data/sample.out", 128);

    Vector rx = func.recover();
    cout << "recover test: " << endl;
    for(int i = 0; i < rx.size(); i++)
        cout << rx[i] << " ";
    cout << endl << endl;

    TriangularInterpolator refunc;
    refunc.load("data/sample.model");
    rx = refunc.recover();
    cout << "load and recover test: " << endl;
    for(int i = 0; i < rx.size(); i++)
        cout << rx[i] << " ";
    cout << endl;
    return 0;
}