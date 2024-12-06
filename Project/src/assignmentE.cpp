#include "BSpline.h"
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

const double sq3 = sqrt(3.0);
const double pi = M_PI;
const double eps = 1e-16;

class Heart{
public:
    Point operator () (const double &t) const{
        return (Point){ sq3*cos(t), 2.0/3.0 * ( sq3*sin(t) + sqrt(sq3*fabs(cos(t))) ) };
    }

    double length(const double &l, const double &r) const{
        double mid = (l+r)/2;
        double A = dist( (*this)(l), (*this)(r) );
        double B = dist( (*this)(l), (*this)(mid) ) + dist( (*this)(mid), (*this)(r) );
        return (fabs(A-B) < eps) ? B : (length(l,mid) + length(mid,r));
    }
} cur;

const double L = -1.5*pi;
const double R = 0.5*pi;
const double len = cur.length(L,R);

int main(int argc, char * argv[]){
    const int n = atoi(argv[2]);
    const double unitlen = len / n; //等距取点时，相邻两点之间曲线段的长度
    cerr << "unit length: " << unitlen << endl;
    vector<double> knots;
    knots.push_back(L);
    double t = L;
    for(int i = 1; i < n; i++){ //二分法计算下一个等距点对应的参数
        double l = t, r = R;
        while( fabs(r-l) > eps*20 ){
            double mid = (l+r)/2;
            if( cur.length(t,mid) < unitlen ) l = mid;
            else r = mid;
        }
        knots.push_back(l);
        t = l;
        cerr << i << " : " << t << endl;
    }
    knots.push_back(R);

    vector<double> x, y;
    for(auto & tv : knots){
        Point tmp = cur(tv);
        x.push_back(tmp.x);
        y.push_back(tmp.y);
    }

    cout << "xp = [ ";
    for(auto & xv : x) cout << xv << " ";
    cout << "];" << endl;
    cout << "yp = [ ";
    for(auto & yv : y) cout << yv << " ";
    cout << "];" << endl;

    Curve splineCurve(x, y, argv[1]);
    cout << splineCurve << endl;

    return 0;
}