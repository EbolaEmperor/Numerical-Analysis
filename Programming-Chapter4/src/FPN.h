#include <vector>
#include <cstdlib>
#include <cmath>

class FPN_System{
public:
    int L, U, b, p;
    FPN_System(){}
    FPN_System(int L, int U, int b, int p):
        L(L), U(U), b(b), p(p) {}
    FPN_System(const FPN_System &fpns):
        L(fpns.L), U(fpns.U), b(fpns.b), p(fpns.p) {}
};

class FPN{
private:
    FPN_System & fpns;
    std::vector<int> m;
    int e;
    bool fg; //fg=1 meas negative
public:
    FPN(FPN_System & _fpns): fpns(_fpns) {
        m.resize(fpns.p);
    }
    FPN(const FPN & rhs): fpns(rhs.fpns), m(rhs.m), e(rhs.e), fg(rhs.fg){}

    void turn_UFLext(){
        for(auto & x : m) x = 0;
        m[fpns.p-1] = 1;
        e = fpns.L;
        fg = 0;
    }

    void turn_UFL(){
        for(auto & x : m) x = 0;
        m[0] = 1;
        e = fpns.L;
        fg = 0;
    }

    void turn_OFL(){
        for(auto & x : m) x = fpns.b-1;
        e = fpns.U;
        fg = 0;
    }

    void turn_zero(){
        fg = 0;
        for(auto & x : m) x = 0;
        e = fpns.L;
    }

    void turn_neg(){
        fg ^= 1;
    }

    bool absIsUFLext(){
        if(e != fpns.L) return false;
        for(int i = 0; i < fpns.p-1; i++)
            if(m[i]) return false;
        return (m[fpns.p-1] == 1);
    }

    bool absIsUFL(){
        if(e != fpns.L) return false;
        for(int i = 1; i < fpns.p; i++)
            if(m[i]) return false;
        return (m[0] == 1);
    }

    bool absIsOFL(){
        bool ismax = true;
        for(auto x : m)
            if(x != fpns.b-1) ismax = false;
        return (ismax && e == fpns.U);
    }

    bool isUFLext(){
        return fg==0 && absIsUFLext();
    }
    bool isUFL(){
        return fg==0 && absIsUFL();
    }
    bool isOFL(){
        return fg==0 && absIsOFL();
    }
    bool isNegUFLext(){
        return fg==1 && absIsUFLext();
    }
    bool isNegUFL(){
        return fg==1 && absIsUFL();
    }
    bool isNegOFL(){
        return fg==1 && absIsOFL();
    }

    bool iszero(){
        for(auto x : m)
            if(x) return false;
        return true;
    }

    void turn_next_abs(){
        int carry = fpns.p - 1;
        for(; carry >= 0 && m[carry] == fpns.b-1; carry--);
        if(carry == -1){
            for(auto & x : m) x = 0;
            m[0] = 1;
            e++;
        } else {
            m[carry]++;
            for(int i = carry+1; i < fpns.p; i++)
                m[i] = 0;
        }
    }

    void turn_pre_abs(){
            int carry = fpns.p - 1;
        for(; carry && m[carry] == 0; carry--);
        if(carry == 0 && m[0]==1){
            for(auto & x : m) x = 1;
            if(e > fpns.L) e--;
            else m[0] = 0;
        } else {
            m[carry]--;
            for(int i = carry+1; i < fpns.p; i++)
                m[i] = fpns.b - 1;
        }
    }

    bool turn_next(){
        if(isOFL()) return false;
        if(iszero()) turn_UFL();
        else if(fg==0) turn_next_abs();
        else if(isNegUFL()) turn_zero();
        else turn_pre_abs();
        return true;
    }

    bool turn_next_ext(){
        if(isOFL()) return false;
        if(iszero()) turn_UFLext();
        else if(fg==0) turn_next_abs();
        else if(isNegUFLext()) turn_zero();
        else turn_pre_abs();
        return true;
    }

    FPN next(){
        FPN tmp(*this);
        if(!tmp.turn_next()) exit(-1);
        return tmp;
    }

    double to_double(){
        double res = 0, base = 1.0;
        for(auto x : m){
            res += x / base;
            base *= fpns.b;
        }
        if(fg) res = -res;
        return res * pow(fpns.b, e);
    }

    FPN operator - () const{
        FPN tmp(*this);
        tmp.fg ^= 1;
        return tmp;
    }
};