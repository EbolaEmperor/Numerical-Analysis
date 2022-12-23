#include "ppForm.h"

using std::cout;
using std::endl;
const double pi = acos(-1);

// Runge函数
class Runge : public Function{
    double operator () (const double &x) const{
        return 1.0/(1.0+x*x);
    }
    double diff (const double &x) const{
        return -2.0*x/pow(1.0+x*x,2);
    }
    virtual double diff2 (const double &x) const{
        return (6.0*x*x-2.0)/pow(1.0+x*x,3);
    }
} runge;

int main(){
    const int n = 7;
    const double xvalue[] = {-5, -3, -1, 0, 1, 3, 5};
    std::vector<double> x(xvalue, xvalue + n);

    ppForm_cubic ppc_natural(x, runge, "natural");
    cout << "y1 = [";
    for(double x = xvalue[0]; x <= xvalue[n-1]; x += 0.01){
        cout << ppc_natural(x) << " ";
    }
    cout << "];" << endl;

    ppForm_cubic ppc_com(x, runge, "complete");
    cout << "y2 = [";
    for(double x = xvalue[0]; x <= xvalue[n-1]; x += 0.01){
        cout << ppc_com(x) << " ";
    }
    cout << "];" << endl;

    ppForm_cubic ppc_nak(x, runge, "not-a-knot");
    cout << "y3 = [";
    for(double x = xvalue[0]; x <= xvalue[n-1]; x += 0.01){
        cout << ppc_nak(x) << " ";
    }
    cout << "];" << endl;

    ppForm_cubic ppc_sec(x, runge, "second-derivatives-at-end");
    cout << "y4 = [";
    for(double x = xvalue[0]; x <= xvalue[n-1]; x += 0.01){
        cout << ppc_sec(x) << " ";
    }
    cout << "];" << endl;

    ppForm_cubic ppc_peri(x, runge, "periodic");
    cout << "y5 = [";
    for(double x = xvalue[0]; x <= xvalue[n-1]; x += 0.01){
        cout << ppc_peri(x) << " ";
    }
    cout << "];" << endl;

    return 0;
}