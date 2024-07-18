#ifndef __PENDULUM_H__
#define __PENDULUM_H__

#include <functional>
#include <iostream>
#include <string>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <TText.h>
#include <TFile.h>
#include <TKey.h>
#include <map>
#include <ODEpoint.h>
#include <Xvar.h>
using namespace std;

class pendulum {
    public:
        pendulum();
        pendulum(double, const Xvar&);
        pendulum(double, const std::initializer_list<double>&); // pendulum(10, {80, 0})
        // length in meters, angles in degrees, velocity degrees/sec
        pendulum(double length = 10, double theta_0 = 80, double theta_vel_0 = 0);
        ~pendulum() = default;

        // solvers
        const std::vector<ODEpoint> StormerVerletSolver(double Time, double step=1E-4);
        const std::vector<ODEpoint> Euler(double Time,double step=1E-4);
        const std::vector<ODEpoint> EulerCromer(double Time,double step=1E-4);
        const std::vector<ODEpoint> Trapezoidal(double Time,double step=1E-4);
        const std::vector<ODEpoint> RungeKutta4(double Time,double step=1E-4);
        const std::vector<ODEpoint> RungeKutta2(double Time,double step=1E-4);
        // graphics
        void Draw(std::string s, string yaxis = "teta", double ti=0., double tf=10.);





    private:
        double L; // length (m)
        Xvar X0; // initial conditions: angle (rad), angular velocity (rad/s)
        // solutions
        std::map<std::string,  std::vector<ODEpoint> > MS; // key: "verlet", "euler", "trapezoidal"
        // functions associated to dependent variables 1st order ODE's
        std::function<double(ODEpoint)> f[2]; // f_theta(t,theta,omega), f_omega(t,theta,omega)
};

#endif