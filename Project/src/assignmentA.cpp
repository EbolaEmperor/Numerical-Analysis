#include "ppForm.h"
using namespace std;

// Runge函数
class Runge : public Function{
public:
    double operator () (const double &x) const{
        return 1.0/(1.0+25*x*x);
    }
    double diff (const double &x) const{
        return -50.0*x/pow(1.0+25*x*x,2);
    }
} runge;

int main(int argc, char * argv[]){
    const int n = atoi(argv[1]);
    ppForm_cubic pp(n, -1, 1+1e-14, runge, "complete");
    cout << "y1 = [ ";
    for(double x = -1; x <= 1.001; x += 0.01)
        cout << pp(x) << " ";
    cout << "];" << endl;

    double maxerr = 0;
    for(int i = 0; i < n-1; i++){
        double l = -1 + 2.0*i/(n-1);
        double r = -1 + 2.0*(i+1)/(n-1);
        double m = (l+r) / 2;
        double err = fabs( runge(m) - pp(m) );
        if(err > maxerr) maxerr = err;
    }
    cerr << "max error: " << maxerr << endl;
    return 0;
}