#include "BSpline.h"
#include "ppForm.h"

using std::cout;
using std::endl;
const double pi = acos(-1);

// 心脏线
class Cardioid_X : public Function{
    double operator () (const double &x) const{
        return 2*sin(x) - sin(2*x);
    }
    double diff (const double &x) const{
        return 2*cos(x) - 2*cos(2*x);
    }
    virtual double diff2 (const double &x) const{
        return -2*sin(x) + 4*sin(2*x);
    }
} hx;

class Cardioid_Y : public Function{
    double operator () (const double &x) const{
        return 2 * ( cos(x) - cos(x)*cos(x) );
    }
    double diff (const double &x) const{
        return 2 * ( -sin(x) + 2*cos(x)*sin(x) );
    }
    virtual double diff2 (const double &x) const{
        return 2 * ( -cos(x) + 2*cos(x)*cos(x) - 2*sin(x)*sin(x) );
    }
} hy;

int main(int argc, char * argv[]){
    const int n = 9;
    const double l = 0, r = 2*pi;

    BSpline_cubic bs_natural_x(n, l, r, hx, argv[1]);
    BSpline_cubic bs_natural_y(n, l, r, hy, argv[1]);
    cout << "x1 = [";
    for(double t = l; t <= r; t += 0.01){
        cout << bs_natural_x(t) << " ";
    }
    cout << "];" << endl;
    cout << "y1 = [";
    for(double t = l; t <= r; t += 0.01){
        cout << bs_natural_y(t) << " ";
    }
    cout << "];" << endl;

    ppForm_cubic ppc_natural_x(n, l, r, hx, argv[1]);
    ppForm_cubic ppc_natural_y(n, l, r, hy, argv[1]);
    cout << "x2 = [";
    for(double t = l; t <= r; t += 0.01){
        cout << ppc_natural_x(t) << " ";
    }
    cout << "];" << endl;
    cout << "y2 = [";
    for(double t = l; t <= r; t += 0.01){
        cout << ppc_natural_y(t) << " ";
    }
    cout << "];" << endl;

   return 0;
}