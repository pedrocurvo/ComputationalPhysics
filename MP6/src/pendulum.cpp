#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <ODEpoint.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <Xvar.h>
#include <pendulum.h>
using namespace std;

pendulum::pendulum(double lenght, const Xvar& X){
    L = lenght;
    X0 = X;
    f[0] = [&](ODEpoint p){return p.X()[1];};
    f[1] = [&](ODEpoint p){return -sin(p.X()[0]);};
}

pendulum::pendulum(double lenght, const std::initializer_list<double>& xvar){
    L = lenght;
    X0 = Xvar(xvar);
    f[0] = [&](ODEpoint p){return p.X()[1];};
    f[1] = [&](ODEpoint p){return -sin(p.X()[0]);};
}

pendulum::pendulum(double length, double theta_0, double theta_vel_0){
    L = length;
    X0 = Xvar({theta_0, theta_vel_0});
    f[0] = [&](ODEpoint p){return p.X()[1];};
    f[1] = [&](ODEpoint p){return -sin(p.X()[0]);};
}

const vector<ODEpoint> pendulum::Euler(double Time,double step){
    int iterations = (int)(Time/step);
    vector<ODEpoint> sol(iterations+1);
    sol[0] = ODEpoint(0, X0);
    for(int i = 1; i <= iterations; i++){
        double omega = sol[i-1].X()[1] + step*f[1](sol[i-1]);
        double teta = sol[i-1].X()[0] + step*f[0](sol[i-1]);
        sol[i] = ODEpoint(i*step, Xvar({teta , omega}));
    }
    MS["euler"] = sol;
    return sol;
}

const vector<ODEpoint> pendulum::EulerCromer(double Time,double step){
    int iterations = (int)(Time/step);
    vector<ODEpoint> sol(iterations+1);
    sol[0] = ODEpoint(0, X0);
    for(int i = 1; i <= iterations; i++){
        double omega = sol[i-1].X()[1] + step*f[1](sol[i-1]);
        double teta = sol[i-1].X()[0] + step * omega;
        sol[i] = ODEpoint(i*step, Xvar({teta , omega }));
    }
    MS["eulercromer"] = sol;
    return sol;
}

const vector<ODEpoint> pendulum::StormerVerletSolver(double Time,double step){
    // int iterations = (int)(Time/step);
    // vector<ODEpoint> sol(iterations+1);
    // for(int i = 0; i < sol.size(); i++){
    //     sol[i] = ODEpoint(i*step, Xvar({0, 0}));
    // }
    // sol[0] = ODEpoint(0, X0);
    // sol[1].X()[0] = X0[0] - step * X0[0] + pow(step, 2) * f[0](ODEpoint(0, X0)) / 2;
    // for(int i = 1; i < iterations; i++){
    //     double teta_i1 = 2 * sol[i].X()[1] - sol[i-1].X()[1] + step*step*f[1](sol[i-1]);
    //     double omega_i = teta_i1 - sol[i-1].X()[1] / 2 / step;
    //     sol[i + 1].X()[0] = teta_i1;
    //     sol[i].X()[1] = omega_i;
    // }
    // MS.insert({"verlet", sol});
    // return sol;
    //Salvador's method
    int iterations = int(Time/step); // numero de espços do vector de ODEpoints
    /*
    vector com ODEpoints: o primeiro espaço do ODEpoint contém o valor de tempo e o seundo espaço 
    contém um Xvar. Cada Xvar contém duas coordeandas: um theta e um omega. Logo, podemos aceder a
    estas variáveis da seguinte forma:
    theta: vec[i].X()[0]
    omega: vec[i].X()[1]
    time: vec[i].T()
    */
    vector <ODEpoint> vec(iterations + 1);
    vec[0] = ODEpoint(0, X0); // o primeiro ODEpoint do vetor contém as condições iniciais
    double theta_n_plus1, theta_n_minus1;

    for (int i = 1; i <= iterations; i++) {
        /*
        Vamos encger cada netrada do vetor com um ODEpoint que tem o tempo correspondente e um
        Xvar que tem um vetor com duas cooordenadas vazias: uma onde iremos colocar o theta e outra
        onde iremos colocar o omega:
        */
        vec[i] = ODEpoint(i * step, Xvar(2));
    }
    /* 
    Vamos preencher o theta do vec[1] porque no loop for seguinte nunca o iremos calcular e
    precisamos do seu valor: 
    */
    vec[1].X()[0] = X0[0] - step*X0[1] + (pow(step, 2))/2 * f[1](vec[0]);
    /*
    Este loop for vai, em cada iteração,
    */
    for (int i = 1; i <= iterations; i++) {
        theta_n_plus1 = 2*vec[i].X()[0] - vec[i - 1].X()[0] + pow(step, 2)*f[1](vec[i]);
        theta_n_minus1 = vec[i - 1].X()[0];
        if (i != iterations) {
            // só colocamos um theta n+1 até à penúltima iteração:
            // theta n plus one:
            vec[i + 1].X()[0] = theta_n_plus1;
        }
        // omega n:
        vec[i].X()[1] = (theta_n_plus1 - theta_n_minus1) / (2*step);
    }
    MS["verlet"] = vec;
    return vec; 
}

const vector<ODEpoint> pendulum::Trapezoidal(double Time,double step){
    int iterations = (int)(Time/step);
    vector<ODEpoint> sol(iterations+1);
    sol[0] = ODEpoint(0, X0);
    for(int i = 1; i <= iterations; i++){
        double omega = sol[i-1].X()[1] + step*f[1](sol[i-1]);
        double teta = sol[i-1].X()[0] + step * (omega + sol[i - 1].X()[1]) / 2;
        sol[i] = ODEpoint(i*step, Xvar({teta , omega }));
    }
    MS["trapezoidal"] =  sol;
    return sol;
};

const std::vector<ODEpoint> pendulum::RungeKutta2(double Time,double step){

	std::vector<ODEpoint> sols;
	sols.push_back(ODEpoint(0,X0));

	double theta = X0.X()[0];
	double theta_ponto = X0.X()[1];
	double t = 0;

	while (t<=Time){
		double meio_t = t + step/2;

		double meio_theta;
		double meio_theta_ponto;

		meio_theta = theta + step/2 * f[0](sols[sols.size()-1]);
		meio_theta_ponto = theta_ponto + step/2 * f[1](sols[sols.size()-1]);

		ODEpoint meio_ponto = ODEpoint(meio_t,Xvar({meio_theta,meio_theta_ponto}));
		t += step;
		theta = theta + step * f[0](meio_ponto);
		theta_ponto = theta_ponto + step * f[1](meio_ponto);

		ODEpoint novo_ponto = ODEpoint(t,Xvar({theta,theta_ponto}));
		sols.push_back(novo_ponto);
	}
	MS["RK2"] = sols;
	return sols;
}

const std::vector<ODEpoint> pendulum::RungeKutta4(double Time, double step) {
    int iterations = int(Time/step);
    vector <ODEpoint> vec(iterations + 1);
    vec[0] = ODEpoint(0, X0);
    // x1 is theta and x2 is omega:
    double k11, k12, k21, k22, k31, k32, k41, k42, x1n, x2n, x1n_plus1, x2n_plus1;
    for (int i = 1; i <= iterations; ++i) {
        x1n = vec[i - 1].X()[0];
        x2n = vec[i - 1].X()[1];
        k11 = f[0](ODEpoint(i*step, Xvar({x1n, x2n})));
        k12 = f[1](ODEpoint(i*step, Xvar({x1n, x2n})));
        k21 = f[0](ODEpoint(i*step + step/2, Xvar({x1n + step/2 * k11, x2n + step/2 * k12})));
        k22 = f[1](ODEpoint(i*step + step/2, Xvar({x1n + step/2 * k11, x2n + step/2 * k12})));
        k31 = f[0](ODEpoint(i*step + step/2, Xvar({x1n + step/2 * k21, x2n + step/2 * k22})));
        k32 = f[1](ODEpoint(i*step + step/2, Xvar({x1n + step/2 * k21, x2n + step/2 * k22})));
        k41 = f[0](ODEpoint(i*step + step, Xvar({x1n + step * k31, x2n + step * k32})));
        k42 = f[1](ODEpoint(i*step + step, Xvar({x1n + step * k31, x2n + step * k32})));

        // fill vector:
        x1n_plus1 = x1n + (step/6)*(k11 + 2*k21 + 2*k31 + k41);
        x2n_plus1 = x2n + (step/6)*(k12 + 2*k22 + 2*k32 + k42);
        vec[i] = ODEpoint(i * step, Xvar({x1n_plus1, x2n_plus1}));
    }
    MS["RK4"] = vec;
    return vec;
}


void pendulum::Draw(string name, string yaxis, double ti, double tf){
    TApplication App("App", nullptr, nullptr);
    TCanvas Canvas("Canvas", "Canvas", 800, 600);
    TMultiGraph* mg = new TMultiGraph();
    auto legend = new TLegend(0.15, 0.15, 0.4, 0.3);
    legend->SetHeader("Legend", "C");


    if(name != "Todos"){
        TGraph* g = new TGraph();
        g->SetLineWidth(4);
        g->SetLineColor(kRed);
        for(auto ODE : MS[name]){
            if(yaxis != "teta"){
                g->AddPoint(ODE.T(), ODE.X()[0]);
            }else{
                g->AddPoint(ODE.T(), ODE.X()[1]);
            }
        };
        legend->AddEntry(g, name.c_str(), "lp");
        mg->Add(g, "lp");
    }else{
        int n = 1;
        for(auto& e : MS){
            TGraph* g2 = new TGraph();
            g2->SetLineWidth(4);
            g2->SetLineColor(n);
            for(auto ODE : e.second){
                if(yaxis != "teta"){
                    g2->AddPoint(ODE.T(), ODE.X()[0]);
                }else{
                    g2->AddPoint(ODE.T(), ODE.X()[1]);
                };
            }
            legend->AddEntry(g2, e.first.c_str(), "lp");
            mg->Add(g2, "pl");
            n++;
        }
    }
    mg->SetTitle("Pendulum");
    mg->GetXaxis()->SetTitle("Time");
    mg->GetYaxis()->SetTitle(yaxis.c_str());
    mg->Draw("a");
    legend->Draw();
    Canvas.Update();
    Canvas.WaitPrimitive();
    gSystem->ProcessEvents();
    // App.Run();
}

