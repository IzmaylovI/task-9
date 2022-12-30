#include <iostream>
#include <iomanip>
#include <fstream>
#include "Runge_Kutta_3.h"
#include <fstream>
#include <clocale>
#include <string>
#include "nlohmann/json.hpp"
#include <windows.h>

using json = nlohmann::json;

std::ofstream data("file.csv");
std::ofstream guide("guide.csv");

//Backup streambuffers of cout 
std::streambuf* stream_buffer_cout = std::cout.rdbuf();
std::streambuf* stream_buffer_cin = std::cin.rdbuf();

//Get the streambuffer of the file
std::streambuf* stream_buffer_data = data.rdbuf();
std::streambuf* stream_buffer_guide = guide.rdbuf();

double maxEn = 0;
double x_mEn = 0;

double maxS = 0;
double x_MS = 0;

double minS = 10000;
double x_mS = 0;

double Xmax = 1000;
int step_control = 0;

double e_down = 1e-7;
double e_up = 5e-4;

double count_mult = 0;
double count_div = 0;

double max_hn = 0;
double x_Mhn = 0;

double min_hn = 1000000;
double x_mhn = 0;

double In = 0;       // The current in the electrical circuit at time t = Xn
double Xn = 0;        // The left boundary of integration (Xn = 0 default)
long long int integr_step = 0;  // Integration number step

void program(const double& L, const double& R, const double& U, const double& I0, const double& H0, const double& epsilon, const double& b) 
{
    //Redirect cout to the data.csv
    std::cout.rdbuf(stream_buffer_data);

    int countr = 0;

    In = I0;       // The current in the electrical circuit at time t = Xn
    Xn = 0;
    double Hn = H0;       // Integration of step
    double _In_1 = 0;     // Numerical value of the current strength at a time t = Xn
    double Ifin = 0;      // 
    double S = 0;         // Estimation of local error
    double mul_s = 0;     // Number of multiplication Hn
    double div_s = 0;     // Number of division Hn
    double err = 0;       // Global error
    double I = I0;        // Analytical value of the current strength at a time t=Xn
    double _In_1prev = -1;

    std::cout.precision(10);

    std::cout << "n, h, t, In, ^In, S, Ifin, I, err, mul_s, div_s" << std::endl;
    std::cout << integr_step<< "," << std::fixed << 0 << "," << Xn << "," << In << "," << _In_1 << "," << S << "," << Ifin << "," << I << "," << err << "," << mul_s << "," << div_s << std::endl;

    integr_step = 1;

    while (Xn <= b && integr_step <= Xmax) 
    {

        if (countr == 600)
            break;
       
        double In_1 = Runge_Kutta_3(Xn, In, Hn, L, R, U);
        
        mul_s = 0;
        
        _In_1 = double_count_half_step(Xn, In, Hn, L, R, U);
        I = reshenie(L, R, U, I0, Xn + Hn);
        
        //S = abs(_In_1 - In_1) / (2 * 2 * 2 - 1);
        S = ((_In_1 - In_1) / (pow(2, 3) - 1));
        err = abs(I - In_1);
        /*
        if (abs(S) > maxS) {
            maxS = abs(S);
            x_MS = Xn;
        }
        if (abs(S) < minS) {
            minS = abs(S);
            x_mS = Xn;
        }

        if (err > maxEn) {
            maxEn = err;
            x_mEn = Xn;
        }
        */
       
        Ifin = In_1 + S;


        if (step_control == 1 || step_control == 2) {
            switch (check(In_1, _In_1, epsilon, e_down, e_up, step_control)) {
            case 0:
                if (abs(S) > maxS) {
                    maxS = abs(S);
                    x_MS = Xn;
                }
                if (abs(S) < minS) {
                    minS = abs(S);
                    x_mS = Xn;
                }

                if (err > maxEn) {
                    maxEn = err;
                    x_mEn = Xn;
                }
                if (Hn > max_hn) {
                    max_hn = Hn;
                    x_Mhn = Xn + Hn;
                }
                if (Hn < min_hn) {
                    min_hn = Hn;
                    x_mhn = Xn + Hn;
                }
                std::cout << integr_step << "," << std::fixed << Hn << "," << std::fixed << Xn + Hn << "," << std::fixed << In_1 << "," << std::fixed << _In_1 << "," << std::fixed << abs(S) << "," << std::fixed << Ifin << "," << std::fixed << I << "," << std::fixed << err << "," << mul_s << "," << div_s << std::endl;
                mul_s = 0;
                div_s = 0;
                Xn = Xn + Hn;
                //In = In_1;
                In = Ifin;
                if (abs(_In_1prev - _In_1) >= Hn) {
                    countr++;
                }
                else {
                    _In_1prev = _In_1;
                    countr = 0;
                }
                
                integr_step++;
                break;
            case 1:
                mul_s++;
                if (abs(S) > maxS) {
                    maxS = abs(S);
                    x_MS = Xn;
                }
                if (abs(S) < minS) {
                    minS = abs(S);
                    x_mS = Xn;
                }

                if (err > maxEn) {
                    maxEn = err;
                    x_mEn = Xn;
                }

                if (abs(_In_1prev - _In_1) >= Hn) {
                    countr++;
                }
                else {
                    _In_1prev = _In_1;
                    countr = 0;
                }

                if (Hn > max_hn) {
                    max_hn = Hn;
                    x_Mhn = Xn + Hn;
                }

                if (Hn < min_hn) {
                    min_hn = Hn;
                    x_mhn = Xn + Hn;
                }

                std::cout << integr_step << "," << std::fixed << Hn << "," << std::fixed << Xn + Hn << "," << std::fixed << In_1 << "," << std::fixed << _In_1 << "," << abs(S) << "," << std::fixed << Ifin << "," << std::fixed << I << "," << std::fixed << err << "," << mul_s << "," << div_s << std::endl;
                Xn = Xn + Hn;
                Hn *= 2;
                //In = In_1;
                count_mult++;
                In = Ifin;
                if (abs(_In_1prev -_In_1) >= Hn) {
                    countr++;
                }
                else {
                    _In_1prev = _In_1;
                    countr = 0;
                }
                integr_step++;
                break;
            case 2:
                Xn = Xn;
                In = In;
                Hn /= 2;
                div_s++;
                count_div++;
                break;
            }
        }
        else {
            /*
            if (abs(S) > maxS) {
                maxS = abs(S);
                x_MS = Xn;
            }
            if (abs(S) < minS) {
                minS = abs(S);
                x_mS = Xn;
            }

            if (err > maxEn) {
                maxEn = err;
                x_mEn = Xn;
            }
            */
            if (Hn > max_hn) {
                max_hn = Hn;
                x_Mhn = Xn + Hn;
            }

            if (Hn < min_hn) {
                min_hn = Hn;
                x_mhn = Xn + Hn;
            }

            Xn = Xn + Hn;
            In = In_1;
            std::cout << integr_step << "," << Hn << ","  << Xn << ","  << In_1 << "," <<  _In_1 << "," << S << ","  << Ifin << ","  << I << "," <<  err << "," << 0 << "," << 0 << std::endl;
        }
    }
}


int main(int argc, char** argv) {
    setlocale(LC_ALL,"rus");
    
    SetConsoleTitleW(L"Команда 4: Измайлов, Митин, Шокуров, Ивлев");

    std::ifstream jsf("k.json");
    json params = json::parse(jsf);

    double L = params["system_param"]["L"];                             // Self - induction coefficient
    double R = params["system_param"]["R"];                             // resistor resistance
    double U = params["system_param"]["U"];                             // battery voltage
    double I0 = params["system_param"]["I0"];                           // Value of the current strength at a time t = 0
    double H0 = params["method_param"]["integration_step"];             // Integration of step
    double epsilon = params["method_param"]["epsilon"];                 // Measurement error
    double b = params["method_param"]["right_border"];                  // The right boundary of Xn
    //e_down = params["method_param"]["e_down"];
    //e_up = params["method_param"]["e_up"];
    std::string control = params["method_param"]["error_control"];
    Xmax = params["method_param"]["max_step"];
    //std::cout << params["step_control"];    // Param for reading control work


//    double epsilon = 1e-4;

    if (control == "No") {
        step_control = 0;
    }
    else if (control == "Up") {
        step_control = 1;
    }
    else {
        step_control = 2;
    }

    //int a = 1;


    //Redirect cout to the data.csv
    std::cout.rdbuf(stream_buffer_cout);
    
    program(L, R, U, I0, H0, epsilon, b);

    std::cout.rdbuf(stream_buffer_guide);
    std::cout.precision(10);
    
    std::cout << "№ варианта задания 5" << std::endl;
    std::cout << "Тип задачи: основная" << std::endl;
    std::cout << "Метод Рунге Кутта порядка 3" << std::endl;
    std::cout << "t0 = ," << 0 << ",I0 = ," << I0 << ",,(Условия задачи Коши)" << std::endl;
    std::cout << "b = ," << b << std::endl;
    std::cout << "h0 = " << std::fixed << H0 << ",Nmax =," << Xmax << std::endl;
    std::cout << "(начальный шаг интегрирования и максимальное число шагов метода)" << std::endl;
    std::cout << "\n";
    std::cout << "epsilon = ," << std::fixed << epsilon << std::endl;
    std::cout << "E_up = ," << epsilon << ",Emin =," << epsilon / pow(2, 4) << std::endl;
    std::cout << "Контроль локальной погрешности сверху и снизу" << std::endl;
    std::cout << "\n";
    std::cout << "Результаты расчета: " << std::endl;
    std::cout << "N =," << integr_step << ",(Количество шагов интегрирования)" << std::endl;
    std::cout << "tитог=," << std::fixed << Xn << ",Iитог=," << In << std::endl;
    std::cout << "\n";
    std::cout << "max|En|=," << std::fixed << maxEn << ",при xn," << std::fixed << x_mEn << std::endl;
    std::cout << "max |S|=," << std::fixed << maxS << ",при xn," << std::fixed << x_MS << std::endl;
    std::cout << "min |S|=," << std::fixed << minS << ",при xn," << std::fixed << x_mS << std::endl;
    std::cout << "\n";
    std::cout << "Всего ум. шага,," << count_div << std::endl;
    std::cout << "Всего ув. шага,," << count_mult << std::endl;
    std::cout << "\n";
    std::cout << "max hn=," << std::fixed << max_hn << ",при x(n+1)," << std::fixed << x_Mhn << std::endl;
    std::cout << "min hn=," << std::fixed << min_hn << ",при x(n+1)," << std::fixed << x_mhn << std::endl;
    //system("ConsoleApplication1.exe");
    system("python graph_copy.py");

    

    return 0;
}