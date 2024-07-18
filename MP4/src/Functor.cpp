#include <iostream>
#include <Functor.h>
using namespace std;

//////////////////////////////
TCanvas* Functor::c = nullptr;

////////////////////////////
//
double Functor::operator()(double x){
    return 1;
}

////////////////////////////
//
void Functor::Draw(double xi, double xf, int num, std::string xtitle, std::string ytitle){
    TGraph* graph = new TGraph(num);
    graph->SetTitle(name.c_str());
    graph->GetXaxis()->SetTitle(xtitle.c_str());
    graph->GetYaxis()->SetTitle(ytitle.c_str());
    double dx = (xf - xi) / (num - 1) ;
    for(int i=0; i<num; i++){
        double x = xi + i * dx;
        graph->SetPoint(i, x, (*this)(x));
       //graph->SetPoint(i, x, operator()(x));

    }
    if (!c) c = new TCanvas("c", "Canvas", 1200, 800);
    graph->Draw("APL");
    graph->SetMarkerStyle(24);
    c->Update();
    c->WaitPrimitive();
    gSystem->ProcessEvents();
}
