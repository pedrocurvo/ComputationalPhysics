#ifndef __INTEGDERIV_H__
#define __INTEGDERIV_H__
#include "MyFunction.h"
using namespace std;

class IntegDeriv {
    public:
        IntegDeriv(Functor&);
        ~IntegDeriv() = default;

        // test function
        void Dump(double, double);
        // integration methods
        void TrapezoidalRule(double xi, double xf, double& Integral, double& Error, bool adaptative = false);
        void simpsonRule(double xi, double xf, double& Integral, double& Error);
        void MonteCarloNormal(double xi, double xf, double& Integral, double& Error);
        void MonteCarloVonNeumann(double xi, double xf, double& Integral, double& Error, int iterations = 1e8);


        // derivative methods
        double SecondDerivative(double x, double h = 1e-5);
        double FirstDerivative(double x0, double h = 1e-5);
        double FourthDerivative(double x, double h = 1e-4);
        double ThirdDerivative(double x, double h = 1e-5);



    private:
        Functor& F;
};

#endif