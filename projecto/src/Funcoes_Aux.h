#ifndef __FUNCOES_AUX_H__
#define __FUNCOES_AUX_H__
#include <iostream>
#include <vector>
#include "TCanvas.h"
using namespace std;

double f(double , double , double , double , double , double );
void SetTemperature(double& , int);
double GetInitialT(double, double, std::vector <std::pair<double,double>>);
void SetParam(vector<double>& vec, vector<double>& vec2);
void Randomize(vector<double>& vec, vector<double>& diff);
vector<double> RandomSystem();
void Calc_Draw_Period(double A,double B,double C,double D, double E,double totaltime,double step, TCanvas* Canvas);



#endif