#include <Functor.h>
#include <MyFunction.h>
#include <iostream>
#include <TApplication.h>
#include <IntegDeriv.h>
using namespace std;

void f(Functor& a){
    TApplication app("app", nullptr, nullptr);
}

int main(){
    MyFunction a("Functor", [&](double x){return abs(sin(x));});
    //TApplication app("app", nullptr, nullptr);
    f(a);
    cout << a.Functor::operator()(3) << endl;
    cout << a.operator()(3) << endl;
    cout << a(3) << endl;
    // a.Draw(0., 10, 10000);
    IntegDeriv f(a);
    double res_2 = f.FirstDerivative(0);
    cout << "First : " << res_2 << endl;
    double res = f.SecondDerivative(0);
    cout << "Second : " << res << endl;
    double res_3 = f.FourthDerivative(0);
    cout << "Fourth : " << res_3 << endl; 
    double Integral = 0;
    double Error = 1e-3;
    f.TrapezoidalRule(0, 1, Integral, Error, true);
    // f.simpsonRule(0, M_PI, Integral, Error);
    //f.MonteCarloNormal(0, M_PI, Integral, Error);
    f.MonteCarloVonNeumann(0, M_PI, Integral, Error);
    cout << "Integral" << endl;
    cout << Integral << endl;
    cout << "Error" << endl;
    cout << Error << endl;
    return 0;
}