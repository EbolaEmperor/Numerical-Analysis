#include "BSpline.h"

using std::cout;
using std::endl;
const double pi = acos(-1);

// Runge函数
class Runge : public Function{
public:
    double operator () (const double &x) const{
        return 1.0/(1.0+x*x);
    }
    double diff (const double &x) const{
        return -2.0*x/pow(1.0+x*x,2);
    }
    double diff2 (const double &x) const{
        return (6.0*x*x-2.0)/pow(1.0+x*x,3);
    }
} runge;

int main(){
    BSpline_quadratic bs_quad(11, -5, 5, runge);
    BSpline_cubic bs_cubic(11, -5, 5, runge, "complete");

    const double ts[] = {-3.5, -3, -0.5, 0, 0.5, 3, 3.5};
    for(int i = 0; i < 7; i++)
        cout << fabs( bs_quad(ts[i]) - runge(ts[i]) ) << "  ";
    cout << endl;
    for(int i = 0; i < 7; i++)
        cout << fabs( bs_cubic(ts[i]) - runge(ts[i]) ) << "  ";
    cout << endl;

    return 0;
}