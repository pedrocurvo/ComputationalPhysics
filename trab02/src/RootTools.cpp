#include "RootTools.h"
using namespace std;

void TDrawVector(double x, double y, double z)
{
    TGraph2D *G_vector = new TGraph2D();
    G_vector->AddPoint(0, 0, 0);
    G_vector->AddPoint(x, y, z);
    G_vector->SetMarkerStyle(22);
    G_vector->SetMarkerColor(kRed);
    G_vector->SetMarkerSize(1.);

    G_vector->SetLineColor(kRed);
    G_vector->Draw("same line p");
}

void CanvasUpdater(TCanvas* c, string name){
    c->Update();
    c->SaveAs(name.c_str());
    c->WaitPrimitive();
    gSystem->ProcessEvents();
}