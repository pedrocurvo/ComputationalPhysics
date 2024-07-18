#include "Funcoes_Aux.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "IntegDeriv.h"
#include <TRandom.h>

using namespace std;

double f(double t, double a, double b, double c, double d, double e){
    double expo = a * exp(-b*t + c);
    double cosseno = cos(d*t - e);
    return expo * cosseno;
}

void SetTemperature(double& T, int k){
    T /= log(k) ;
}

double GetInitialT(double S, double X0, vector <pair<double,double>> t_teta_m){
    vector<pair<double, double>> sistemas;
    for(int i = 0; i<sistemas.size(); ++i){
        double a = gRandom -> Rndm();
        double b = gRandom -> Rndm();
        double c = gRandom -> Rndm();
        double d = gRandom -> Rndm();
        double e = gRandom -> Rndm();
        double a_2 = gRandom -> Rndm();
        double b_2 = gRandom -> Rndm();
        double c_2 = gRandom -> Rndm();
        double d_2 = gRandom -> Rndm();
        double e_2 = gRandom -> Rndm();
        double sum_1 = 0.0;
        double sum_2 = 0.0;
        for(int j = 0; j<S; ++j){
            sum_1 += t_teta_m[i].second - f(t_teta_m[i].first,a,b,c,d,e);
            sum_2 += t_teta_m[i].second - f(t_teta_m[i].first,a_2,b_2,c_2,d_2,e_2);
        }
        sistemas.push_back(make_pair(sum_1, sum_2));
    }
    double T1=0;
    for(auto el : sistemas){
        T1 -= abs(el.first-el.second);
    }
    T1 /= (S*log(X0));
    double T_final = 0;
    while(true){
        double num = 0;
        double den = 0;
        for(auto el : sistemas){
            double max = (el.first > el.second)? el.first : el.second;
            double min = (el.first < el.second)? el.first : el.second;
            num += exp(-max/T1);
            den += exp(-min/T1);
        }
        double x_chapeu = num/den;
        if(abs(x_chapeu-X0) < 1e-3){
            T_final = T1;
            break;
        }
        T1 *= pow(log(x_chapeu)/log(X0), 1/2); // p=2
    }
    return T_final;
}

void SetParam(vector<double>& vec, vector<double>& vec2){
    for(int i = 0; i < vec.size(); i++){
        vec[i] = vec2[i];
    };
}

void Randomize(vector<double>& vec, vector<double>& diff){
    for(int i = 0; i < vec.size(); i++){
        vec[i] -= diff[i] * gRandom -> Rndm();
    };
}

vector<double> RandomSystem(){
    vector<double> vec;
    double random = gRandom -> Rndm();
        vec.push_back(gRandom->Gaus(2, 0.1));
        vec.push_back(gRandom->Gaus(0.005, 2));
        vec.push_back(gRandom->Gaus(0, 2));
        vec.push_back(gRandom->Gaus(2, 2));
        vec.push_back(gRandom->Gaus(0, 2));

    return vec;
}

void Calc_Draw_Period(double A,double B,double C,double D, double E,double totaltime,double step, TCanvas* Canvas){
    vector < pair <double,double> > periodos;
    auto teta = [&](double tempo){
        return A*exp(-B * tempo+ C)*cos(D*tempo - E);
    };
    MyFunction tetafunc("Functor",teta);
    IntegDeriv funcao(tetafunc);
    double t = 0;
    while (true){
        //condicao para aplicar formula -> velocidade igual a zero
        //cout << "antes" << endl;
        if (abs(funcao.FirstDerivative(t)) < 1e-1){
            double ke = sin(teta(t) * 0.5);
            auto f = [&](double u){
                return 1/(sqrt(1 - ke * ke * sin(u)*sin(u)));
            };
        //calcular integral
        MyFunction a("Functor",f);
        IntegDeriv fun(a);
        //integrar f
        double Error = 1e-4;
        double Integral = 0;
        fun.TrapezoidalRule(0, 2*M_PI, Integral, Error);

        //calcular periodo
        double T = sqrt(2.6 / 9.8) * Integral;
            periodos.push_back(make_pair(t,T));
        }

        t+=step;
        if (t>=totaltime){
            break;
        }
    }
    
    TGraph* g = new TGraph();
    g->SetTitle("Period");
    g->GetXaxis()->SetTitle("Time(s)");
    g->GetYaxis()->SetTitle("Period");
    g->SetLineWidth(4);
    g->SetLineColor(kOrange+7);
    g->SetMarkerColor(kOrange+7);
    for(pair <double,double> par : periodos){
        g->AddPoint(par.first, par.second);
    }
    g->GetXaxis()->SetRangeUser(0, totaltime);
    g->Draw("APL");
    Canvas->Update();
    Canvas->SaveAs("periodo_medida.png");
    Canvas->WaitPrimitive();
    gSystem->ProcessEvents();
    Canvas->Clear();
}



















