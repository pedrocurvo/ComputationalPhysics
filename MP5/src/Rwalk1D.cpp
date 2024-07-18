#include "Rwalk1D.h"
#include <vector>
#include <map>
#include <TGraph.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TAxis.h>
#include <TSystem.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <TRandom.h>
using namespace std;

Rwalk1D::Rwalk1D(int n, double x, double pl, double pr, double dT, double dX)
{
    x0 = x;
    N = n;
    pL = pl;
    pR = pr;
    dt = dT;
    dx = dX;
    for (int i = 0; i < N; i++){
        vector<double> vec;
        vec.push_back(x0);
        mT.insert({i, vec});
    }
}

const std::vector<double>& Rwalk1D::GetTrajectory(int n){
    return mT[n];
}

double Rwalk1D::GetTimeStep(){
    return dt;
}

double Rwalk1D::GetSpaceStep(){
    return dx;
}
double Rwalk1D::GetNumber(){
    return N;
}

void Rwalk1D::Run(int nsteps){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < nsteps; j++){
            double n = gRandom->Rndm();
            if(n >= pL){
                mT[i].push_back(mT[i].back() + dx);
            }else{
                mT[i].push_back(mT[i].back() - dx);
            }
        }
    }
}

void Rwalk1D::Draw(TCanvas* c){
    int n  = 0;
    double max = 0;
    double min = 0;
    for(int i = 0; i < N; i++){
        for (int j = 0; j < mT[i].size(); j++){
            if(mT[i][j] > max){
                max = mT[i][j];
            }
            if(mT[i][j] < min){
                min = mT[i][j];
            }
        }
    }
    
    for (int i = 0; i < N; i++){
        n++;
        TGraph* g = new TGraph;
        g->SetTitle("Random Walk");
        g->GetXaxis()->SetTitle("time");
        g->GetYaxis()->SetTitle("x");
        g->SetLineColor(n);
        g->GetYaxis()->SetRange(0, 10000);
        g->SetLineWidth(4);
        for (int j = 0; j < mT[i].size(); j++){
            g->AddPoint(j*dt, mT[i][j]);
        }
        if(i == 0){
            g->Draw();
            g->GetYaxis()->SetRangeUser(min, max);
            c->Update();
        }else{
            g->Draw("same");
            g->GetYaxis()->SetRangeUser(min, max);
            c->Update();
        }
    }
    c->WaitPrimitive();
    gSystem->ProcessEvents();
    
}