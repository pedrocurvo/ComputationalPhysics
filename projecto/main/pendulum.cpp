/*****************************************/
/*                                       */
/*            21-24/04/2022              */
/*              'Projeto'                */
/*             Group: D01                */
/*                                       */
/*****************************************/
#include <Xvar.h>
#include <ODEpoint.h>
#include "ODEsolver.h"
#include <RootTools.h>
#include <pendulum_p.h>
#include <iostream>
using namespace std;

int main(){

    double k = 0.01; // coeficiente de amortecimento
    double L = 2.60; // lenght of pendulum
    double m = 0.5; // mass of pendulum
    double g = 9.81; // acceleration of gravity
    double A = k * sqrt(L / g);
    //Vector with ODEs
    vector<std::function<double(ODEpoint)>> functions(2);
    //Placing funtions (eq diferenciais de primeira ordem obtidas a partir da eq do movimento)
    auto f0 = [&](ODEpoint p){return p.X()[1];};
    auto f1 = [&](ODEpoint p){return - A * p.X()[1] - sin(p.X()[0]);};
    functions[0] = f0;
    functions[1] = f1;

    //Initial Conditions teta(0) = 70 graus; teta_ponto(0) = 0 rad/s; coeficiente de amortecimento = 0.01s^-1
    double angle = 70; //graus
    double angle_rad = angle * M_PI / 180; // change to rad
    // the second 0 is the initial angular velocity
    ODEpoint p(0, {angle_rad, 0});
    ODESolver pend_ode(functions, p);

    //Total Time = 200s
    double time = 200;
    //Step = 0.01s
    double step = 0.01;
    //Solving (variables are multiplied by sqrt(g / L) to get the correct time step for non dimensional solver)
    vector<ODEpoint> RK4 = pend_ode.RungeKutta4(time * sqrt(g / L), step * sqrt(g / L));
    
    //Creating Pendulum
    pendulum pendulum(RK4, L, m, k);

    //Graphs
    TApplication App("App", nullptr, nullptr);
    auto Canvas = new TCanvas("Canvas", "Canvas", 1200, 900);

    //Graph teta(t)
    pendulum.DrawTeta(Canvas, time);
    //Graph teta_ponto(t).vs.teta(t)
    pendulum.DrawTetaVSTetaPoint(Canvas, time);
    pendulum.DrawTetaANDTetaPoint(Canvas, time);
    //Total Energy of Pendulum
    pendulum.DrawEnergy(Canvas, time);
    //Expression of Tension Force magnitude
    pendulum.DrawTension(Canvas, 50);
    //Sobreposition of 2 periods over time (calculated with amplitude and time)
    pendulum.DrawPeriod(Canvas, time);

    //Dissipated Energy of Pendulum over 50s with k1 = 0.5kf;
    double res = pendulum.DissipatedEnergy(50);
    cout << "Energy Dissipated from Fa after 50s : " << res << " J" << endl;
    return 0;
}








