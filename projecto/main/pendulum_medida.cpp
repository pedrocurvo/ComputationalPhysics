#include "ODEsolver.h"
#include "pendulum_p.h"
#include <iostream>
#include <TRandom.h>
#include <random>
#include <cmath>
#include "Funcoes_Aux.h"
using namespace std;

int main() {

	double k = 0.01;
    double L = 2.6;
    double m = 0.5;
    double g = 9.8;

    //Vector with ODEs
    vector<std::function<double(ODEpoint)>> functions;
    //Placing funtions
    double A = sqrt(L / g) * k;
    auto f0 = [&](ODEpoint p){return p.X()[1];};
    auto f1 = [&](ODEpoint p){return - A * p.X()[1] - sin(p.X()[0]);};
    functions.push_back(f0);
    functions.push_back(f1);

    //Initial Conditions teta(0) = 70 graus; teta_ponto(0) = 0 rad/s; coeficiente de amortecimento = 0.01s^-1
    double angle = 70; //graus
    double angle_rad = angle * M_PI / 180;
    ODEpoint p(0, {angle_rad, 0});
    ODESolver pend_ode(functions, p);

    //Total Time = 200s
    double time = 200;
    //Step = 0.01s
    double step = 0.05;
    //Solving
    // multiply by sqrt(g / L) to get the correct time step for non dimensional solver 
    vector<ODEpoint> RK4 = pend_ode.RungeKutta4(time * sqrt(g/L), step * sqrt(g/L));

    //Creating Pendulum
    pendulum pendulum(RK4, L, m, k);

    //get dimensional ode points
    vector<ODEpoint> odes = pendulum.GetPointsDimension();

    //create vector of pairs <t,teta>
    vector <pair<double,double>> t_teta; 
    for (ODEpoint P : odes){
    	t_teta.push_back(make_pair(P.T(),P.X()[0]));
    }

    //get teta_m for first 50s -> vector <t,teta_m>
    vector <pair<double,double>> t_teta_m;
    int i = 0;
    double temp = t_teta[i].first;
    while (temp < 50.){
    	double G = gRandom->Gaus(0,0.05 * t_teta[i].second);
    	double teta_m = t_teta[i].second + G;
    	t_teta_m.push_back(make_pair(t_teta[i].first,teta_m));
    	cout << "tempo: " << t_teta[i].first << " | " << "teta : " << t_teta[i].second  << " | " << "teta_m: " << teta_m  << endl;
    	++i;
    	temp = t_teta[i+1].first;
    }

    //calcular dispersao relativa para cada t
    vector <pair<double,double>> disp;
    for (int i = 0; i < t_teta_m.size(); i++){
    	double d = (t_teta_m[i].second - t_teta[i].second)/t_teta[i].second;
    	disp.push_back(make_pair(t_teta[i].first , d));
    } 

    //Graphs
    TApplication App("App", nullptr, nullptr);
    auto Canvas = new TCanvas("Canvas", "Canvas", 1200, 900);

    //Graph teta(t)
    pendulum.DrawTeta(Canvas, 50,"teta_medida.png");

    //Graph Teta_m(t)
    TGraph* gteta_m = new TGraph();
    gteta_m->SetTitle("teta_m(t)");
    gteta_m->GetXaxis()->SetTitle("tempo(s)");
    gteta_m->GetYaxis()->SetTitle("teta_m");
    gteta_m->SetLineWidth(4);
    gteta_m->SetLineColor(kAzure + 5);
    gteta_m->SetMarkerColor(kAzure + 5);
    for(pair <double,double> par : t_teta_m){
        gteta_m->AddPoint(par.first,par.second);
    }
    gteta_m->Draw("APL");
    Canvas->Update();
    Canvas->SaveAs("Teta_m(t).png");
    Canvas->WaitPrimitive();
    gSystem->ProcessEvents();
    Canvas->Clear();

    //draw dispersao relativa(t)
    TGraph* graph = new TGraph();
    graph->SetTitle("dispersao_relativa(t)");
    graph->GetXaxis()->SetTitle("tempo(s)");
    graph->GetYaxis()->SetTitle("dispersao_relativa");
    graph->SetLineWidth(1);
    graph->SetMarkerSize(1);
    graph->SetMarkerStyle(43);
    graph->SetMarkerColor(kGreen + 2);
    graph->SetLineColor(kGreen + 2);
    for(pair <double,double> par : disp){
        graph->AddPoint(par.first,par.second);
    }
    graph->Draw("APL");
    Canvas->Update();
    Canvas->SaveAs("Teta.png");
    Canvas->WaitPrimitive();
    gSystem->ProcessEvents();
    Canvas->Clear();

    //NEXT EXERCISE
    // Alínea (b)
    double X0 =0.6;
    double S = 100/X0; // S tem de ser grande
    // double T0 = GetInitialT(S,X0,t_teta_m);
    double T0 = 2;
    // Parametros iniciais
    // a b c d e 
    vector <double> param = RandomSystem();
    int contador = 2;
    auto const hes = std::random_device{}();
    gRandom->SetSeed(hes);
    int it = 1000;
    double erro = 0;
    // Heavy algorithm, it will take a few seconds to run and draw the graph

    while(it < 1e5){
        vector<double> param2 = RandomSystem(); // Shake system with random values 
        double E1 = 0;
        double E2 = 0;
        erro = 0;
        for(int i = 0; i<t_teta_m.size(); ++i){
            E1 += pow(t_teta_m[i].second - f(t_teta_m[i].first, param[0],param[1],param[2],param[3],param[4]), 2);
            erro += t_teta_m[i].second - f(t_teta_m[i].first, param[0],param[1],param[2],param[3],param[4]);
            E2 += pow(t_teta_m[i].second - f(t_teta_m[i].first, param2[0],param2[1],param2[2],param2[3],param2[4]), 2);
        }
        double diferenca = E2 - E1; // See Error
        if (diferenca < 0){
            SetTemperature(T0,contador); // Change T
            SetParam(param, param2); // Change System
        }else{
            double probability = exp(- diferenca / T0);
            double random = gRandom -> Rndm();
            if(probability > random){
                SetTemperature(T0,contador); // Change T
                SetParam(param, param2); // Change System
            }
        }
        contador++; 
        it++;
    }
    // Amplitude
    for(int i = 0; i < 1e4; i++){
        double A_initial = param[0];
        double A2_initial = gRandom->Gaus(param[0], 0.1);
        double E1 = 0;
        double E2 = 0;
        for(int i = 0; i<t_teta_m.size(); ++i){
            E1 += pow(t_teta_m[i].second - f(t_teta_m[i].first, A_initial,param[1],param[2],param[3],param[4]), 2);
            E2 += pow(t_teta_m[i].second - f(t_teta_m[i].first, A2_initial,param[1],param[2],param[3],param[4]), 2);
        }
        double diferenca = E2 - E1; // See Error
        if (diferenca < 0){
            param[0] = A2_initial;
        }else{
            double probability = exp(- diferenca / T0);
            double random = gRandom -> Rndm();
            if(probability > random){
                param[0] = A2_initial;
            }
        }
    }
    // Great number of amplitudes error might be zero
    cout << endl << "Erro (soma dos quadrados minimos / N_points) : " << abs(erro) / t_teta_m.size() << endl << endl;
    cout << endl << "T0" << T0 << endl;
    cout << "Parametros " << param[0] << " " << param[1] << " " << param[2] << " " << param[3] << " " << param[4] << endl;
    // Be aware that the Graph might get strange, this alghoritm is meant to be used several times until they overlap
    // After various iterations, the graph will be more accurate, at least 3 runtimes
    pendulum.DrawTetaTF1(Canvas, 50, param[0],param[1],param[2],param[3],param[4]);
    Calc_Draw_Period(param[0],param[1],param[2],param[3],param[4], 200 ,0.05 ,Canvas);
    return 0;
}

















