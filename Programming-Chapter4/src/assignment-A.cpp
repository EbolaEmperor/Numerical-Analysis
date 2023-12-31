#include <iostream>
#include <cmath>

double f(const double &x){
    return pow(x,8) - 8*pow(x,7) + 28*pow(x,6) - 56*pow(x,5) + 70*pow(x,4) - 56*pow(x,3) + 28*x*x - 8*x + 1;
}

double g(const double &x){
    return ( ( ( ( ( ( ( x - 8 )*x + 28 )*x - 56 )*x + 70 )*x - 56 )*x + 28 )*x - 8 )*x + 1;
}

double h(const double &x){
    return pow(x-1, 8);
}

const double L = 0.99, R = 1.01;
const double STEP = 0.0002;

int main(){
    for(double x = L; x <= R; x += STEP){
        std::cout.width(15);
        std::cout << f(x);

        std::cout.width(15);
        std::cout << g(x);

        std::cout.width(15);
        std::cout << h(x);
        
        std::cout << ";" << std::endl;
    }
    return 0;
}