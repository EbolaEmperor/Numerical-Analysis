#ifndef __ARNOLDI_INTERP_H__
#define __ARNOLDI_INTERP_H__

#include <cassert>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
using namespace Eigen;

class ArnoldiInterp{
private:
    VectorXd d;
    MatrixXd H;
public:
    ArnoldiInterp(const VectorXd &x,
                  const VectorXd &f){
        int n = x.size() - 1;
        H = MatrixXd::Zero(n + 1, n);
        MatrixXd Q = MatrixXd::Zero(n + 1, n + 1);
        Q.col(0).setConstant(1);
        
        for(int k = 0; k < n; k++){
            VectorXd q = x.cwiseProduct(Q.col(k));
            for(int j = 0; j <= k; j++){
                H(j, k) = Q.col(j).dot(q / n);
                q -= H(j, k) * Q.col(j);
            }
            H(k + 1, k) = q.norm() / sqrt(n);
            Q.col(k + 1) = q / H(k + 1, k);
        }
        
        d = Q.colPivHouseholderQr().solve(f);
    }

    ArnoldiInterp(const std::vector<double> &_x,
                  const std::vector<double> &_f) {
        int n = _x.size() - 1;
        VectorXd x(n + 1), f(n + 1);
        for(int i = 0; i <= n; i++) {
            x(i) = _x[i];
            f(i) = _f[i];
        }
        *this = ArnoldiInterp(x, f);
    }

    VectorXd operator () (const VectorXd &s) const {
        int M = s.size(), n = d.size() - 1;
        MatrixXd W = MatrixXd::Zero(M, n + 1);
        W.col(0).setConstant(1);

        for (int k = 0; k < n; k++) {
            VectorXd w = s.cwiseProduct(W.col(k));
            for (int j = 0; j <= k; j++) {
                w -= H(j, k) * W.col(j);
            }
            W.col(k + 1) = w / H(k + 1, k);
        }

        return W * d;
    }

    std::vector<double> operator () (const std::vector<double> &s) const {
        int M = s.size();
        VectorXd sv(M);
        for(int i = 0; i < M; i++)
            sv(i) = s[i];
        VectorXd fv = operator()(sv);
        std::vector<double> f(M);
        for(int i = 0; i < M; i++)
            f[i] = fv(i);
        return f;
    }

    double operator () (const double &x) const {
        VectorXd s(1);
        s(0) = x;
        return operator()(s)(0);
    }
};


#endif