#include "BSpline.h"
using namespace std;

double RAND(){
    return (double)rand()/RAND_MAX;
}

int main(){
    srand(time(0));
    const int n = 5;
    vector<double> x, y;
    for(int i = 0; i < n; i++){
        x.push_back(RAND()*20);
        y.push_back(RAND()*20);
    }
    x.push_back(x.front());
    y.push_back(y.front());

    cout << "xp = [ ";
    for(auto & v : x) cout << v << " ";
    cout << "];" << endl;
    cout << "yp = [ ";
    for(auto & v : y) cout << v << " ";
    cout << "];" << endl << endl;

    Curve cve1(x, y, "natural");
    cout << cve1 << endl;

    Curve cve2(x, y, "periodic");
    cout << cve2 << endl;

    return 0;
}