#include "Runge_Kutta_3.h"

double Runge_Kutta_3(const double& Xn, const double& In, const double& Hn, const double& L, const double& R, const double& U) {

    double k1 = function(Xn, In, L, R, U);
    double k2 = function(Xn + Hn / 2, In + Hn / 2 * k1, L, R, U);
    double k3 = function(Xn + Hn, In + Hn * (-k1 + 2 * k2), L, R, U);

    double In_1 = In + Hn / 6 * (k1 + 4 * k2 + k3);


    return In_1;
 }

double function(const double& Xn, const double& In, const double& L, const double& R, const double& U) {
    return (U - In * R) / L;
}

double double_count_half_step(double Xn, double In, double Hn, double L, double R, double U) {
    Hn /= 2;

    double _In_1 = Runge_Kutta_3(Xn, In, Hn, L, R, U);
    Xn += Hn;
    _In_1 = Runge_Kutta_3(Xn, _In_1, Hn, L, R, U);

    return _In_1;
}                                                                                                 

int check(double In_1, double _In_1, const double& epsilon, double& e_down, double& e_up, int& step_control) {
    double answer = abs(_In_1 - In_1) / (pow(2, 3) - 1);
    double eps1 = epsilon / (2*2*2*2);
    double eps2 = epsilon;                                                                                                                                                                                                                             
    
    int ans = 0;
    
    if (step_control == 1) {
        if (answer >= eps2) {
            ans = 2;
        }
    }
    else {
        if (answer >= eps2) {
            ans = 2;
        }
        else if (answer <= eps1) {
            ans = 1;
        }
    }

    /*
    if (answer >= eps1 && answer <= eps2) {
        ans = 0;
    }
    else if (answer < eps1) {
        ans = 1;
    }
    else{
        ans = 2;
    }
    */
    /*
    int ans = 0;

    if (step_control == 1) {
        if (answer > e_up) {
            ans = 2;
        }
    }
    else {
        if (answer > e_up) {
            ans = 2;
        }
        else if (answer < e_down) {
            ans = 1;
        }
    }
    */

    return ans;
}

double reshenie(const double& L, const double& R, const double& U, const double& I0, const double& Xn) {
    double st = R * Xn / L;
    double k = exp(st);
    double m = (U - R * I0);
    double mm = (U - R * I0) / exp(st);
    double mmm = U - mm;
    double ans = (double)((U  - (U - R * I0) / exp(st)) / R);
    if (isnan(ans)) {
        std::cout << ",,,,," << R << ',' << Xn << ',' << L << ',' << st << std::endl;
    }
    return ans;
}