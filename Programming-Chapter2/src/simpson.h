#ifndef __SIMPSON_H__
#define __SIMPSON_H__

double simpson(auto f, double a, double b) {
    return (b - a) / 6.0 * (f(a) + 4.0 * f((a + b) / 2.0) + f(b));
}

double adpt_simpson(auto f, double a, double b, double eps, double whole) {
    double c = (a + b) / 2.0;
    double left = simpson(f, a, c);
    double right = simpson(f, c, b);
    if (std::abs(left + right - whole) <= 15.0 * eps)
        return left + right + (left + right - whole) / 15.0;
    return adpt_simpson(f, a, c, eps / 2.0, left) + adpt_simpson(f, c, b, eps / 2.0, right);
}

double adpt_simpson(auto f, double a, double b, double eps = 1e-6) {
    return adpt_simpson(f, a, b, eps, simpson(f, a, b));
}

#endif