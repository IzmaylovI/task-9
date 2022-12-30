#pragma once
#include <cmath>
#include <iostream>


double Runge_Kutta_3(const double& Xn, const double& In, const double& Hn, const double& L, const double& R, const double& U);
double function(const double& Xn, const double& In, const double& L, const double& R, const double& U);

double double_count_half_step(double Xn, double In, double Hn, double L, double R, double U);

double reshenie(const double& L, const double& R, const double& U, const double& I0, const double& Xn);

//int check(double In_1, double _In_1, const double& epsilon);
int check(double In_1, double _In_1, const double& epsilon, double& e_down, double& e_up, int& step_control);