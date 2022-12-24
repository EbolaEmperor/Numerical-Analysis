#include "BSpline.h"

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
    BSpline_quadratic bs_quad(11, -5, 5, runge);
    cout << "y1 = [ ";
    for(double x = -5; x <= 5; x += 0.01)
        cout << bs_quad(x) << " ";
    cout << "];" << endl;

    BSpline_cubic bs_cubic(11, -5, 5, runge, "complete");
    cout << "y2 = [ ";
    for(double x = -5; x <= 5; x += 0.01)
        cout << bs_cubic(x) << " ";
    cout << "];" << endl;

    return 0;
}