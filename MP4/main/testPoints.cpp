#include <DataPoints.h>
#include <Interpolator.h>
using namespace std;

int main(){
    //DataPoints P({{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    //cout << P << endl;
    // P.Draw();

    //Interpolator I({{-2, 5}, {1, 7}, {3, 11}, {7, 34}});
    Interpolator I({{1, 0}, {2, 1}, {3, 0}, {4, 1}, {5, 0}});
    cout << I << endl;
    double var = I.InterpolateLagrange(3.5);
    double var_newton = I.InterpolateNewton(3.5);
    double var_spline = I.InterpolateSpline3(3.5);
    cout << var << endl;
    cout << var_newton << endl;
    cout << var_spline << endl;
    I.Draw("Todos");
    //I.Draw("Todos");
    
    return 0;
}