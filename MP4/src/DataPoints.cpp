#include "DataPoints.h"
#include <iostream>
#include <stdexcept>
using namespace std;

DataPoints::DataPoints(int N, double * x, double * y)
{
    for(int i = 0; i < N; i++)
    {
        P.emplace_back(x[i], y[i]); // se der erro meter makepair
    }
}

DataPoints::DataPoints(const vector<double> & x, const vector<double> & y)
{   
    if(x.size() == y.size())
    {
        for(int i = 0; i < x.size(); i++)
        {
            P.emplace_back(x[i], y[i]);
        }
    }else{
        throw invalid_argument("x and y must have the same size");
    }
}

DataPoints::DataPoints(const vector<pair<double, double>> & v) : P(v) {;}

DataPoints::~DataPoints()
{
    P.clear();
}

//Getters
const vector<pair<double,double> >& DataPoints::GetPoints() const
{
    return P;
}

void DataPoints::GetGraph(TGraph& g)
{
    for(int i = 0; i < P.size(); i++)
    {
        g.AddPoint(P[i].first, P[i].second);
    }
    g.SetMarkerStyle(25);
    g.SetMarkerSize(2.);
    g.SetMarkerColor(kBlue+2);
}

void DataPoints::Draw()
{
    TGraph g;
    GetGraph(g);
    TApplication A("Data Points", nullptr, nullptr);
    if(!c) c = new TCanvas("c", "Data Points", 0, 0, 1200, 900);
    g.Draw("AP");
    c->Update();
    c->WaitPrimitive();
    gSystem->ProcessEvents();
    
}

std::ostream& operator<<(std::ostream& s, const DataPoints& A){
    s << "[";
    for(int i = 0; i < A.P.size(); i++)
    {
        s << "(" << A.P[i].first << ", " << A.P[i].second << ")"<< endl;
    }
    s << "]";
    return s;
}

//TCanvas
TCanvas* DataPoints::c = nullptr;