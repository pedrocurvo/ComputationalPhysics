#include <iostream>
#include <vector>
#include <ODEpoint.h>
#include <Xvar.h>
#include <pendulum.h>
using namespace std;

int main(){

    ODEpoint p(2.6, {70, 0});
    //cout << p << endl;
    pendulum pend(10, 80, 0);
    const vector<ODEpoint>& rungekutta4 = pend.RungeKutta4(200, 0.1);

    
    cout << "Runge Kutta4" << endl << endl;
    for(int i = 0; i < rungekutta4.size(); i++){
        cout << rungekutta4[i] << endl;
    }
    cout << "Graph" << endl << endl;
    pend.Draw("RK4", "teta", 0, 100);
    return 0;
}