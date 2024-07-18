#include <iostream>
#include <functional>
#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <TText.h>
#include <TFile.h>
#include <TKey.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
using namespace std;

int main(){
    double t;
    double x[2];
    // constanța
    double g = 9.8;
    //setting the functions and variables
    function<double(double, double, double)> f[2];
    f[0] = [](double t, double x, double v){return v;};
    f[1] = [g](double t, double x, double v){return g;};

    //setting the initial conditions
    double t0 = 0;
    double tmax = 10;
    x[0] = 0;
    x[1] = 0;
    double step = 1e-1;
    t = t0;
    //vector of Coordinates
    vector<vector<double>> trajectory;
    trajectory.push_back({t, x[0], x[1]});

    while(t<tmax){
        //calculating the next point
        // x(t+1) = x(t) + v(t) * h
        x[0] = x[0] + step * f[0](t, x[0],x[1]);
        x[1] = x[1] + step * f[1](t, x[0],x[1]);
        t += step;
        vector<double> point {t, x[0], x[1]};
        trajectory.emplace_back(point);
    }

    //making the plot
    auto graph = new TGraph();
    for(auto x : trajectory){
        graph->AddPoint(x[0], x[1]);
        cout << "1: " << x[0] << "2 : " << x[1] << endl;
    }
    TApplication a("a", 0, 0);
    TCanvas c("c", "c", 1200, 900);
    graph->Draw();
    c.Update();
    c.WaitPrimitive();

    return 0;
}