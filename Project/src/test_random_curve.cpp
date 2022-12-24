#include "BSpline.h"
using namespace std;

struct Point{
    double x, y;
};

class Curve{
private:
    int n;
    BSpline_cubic fx, fy;

public:
    Curve() : n(0) {}
    Curve(const vector<double> &x, const vector<double> &y, const string & bondary){
        if(x.size() != y.size()) error("Curve :: vector x and y must be the same length!");
        n = x.size();
        vector<double> t;
        for(int i = 0; i < x.size(); i++)
            t.push_back(i);
        fx = BSpline_cubic(t, x, bondary);
        fy = BSpline_cubic(t, y, bondary);        
    }

    Point operator () (const double &t){
        if(t < 0 || t > 1) error("Curve :: Out of range [0,1]!");
        if(n == 0) error("Curve :: Visited operator () at an empty curve!");
        return (Point){fx(t*(n-1)), fy(t*(n-1))};
    }

    friend ostream & operator << (ostream & out, const Curve &rhs){
        static int cnt = 0;
        cnt++;
        out << "x" << cnt << " = [ ";
        for(double i = 0; i <= rhs.n-1; i += 0.01)
            out << rhs.fx(i) << " ";
        out << "];" << endl;
        out << "y" << cnt << " = [ ";
        for(double i = 0; i <= rhs.n-1; i += 0.01)
            out << rhs.fy(i) << " ";
        out << "];" << endl;
        return out;
    }
};

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