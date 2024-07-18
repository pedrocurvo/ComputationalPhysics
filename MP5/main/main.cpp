#include "Rwalk1D.h"
#include "TRandom.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TApplication.h>
using namespace std;

int main() {

    Rwalk1D mov(50, 0, 0.5, 0.5, 1, 1);
    TApplication A("Data Points", nullptr, nullptr);
    mov.Run(500);
    auto c = new TCanvas("c", "Data Points", 0, 0, 1200, 900);
    mov.Draw(c);


    // // ---------------------------
    cout << endl << "Médias:" << endl << endl;
    Rwalk1D mov1(200);
    mov1.Run(10);
    Rwalk1D mov2(200);
    mov2.Run(100);   

    double sum = 0;
    double max1 = 0;
    for(int i = 0; i < mov1.GetNumber(); ++i){
        sum += mov1.GetTrajectory(i).back() - mov1.GetTrajectory(i).front();
        if (mov1.GetTrajectory(i).back() > max1){max1 = mov1.GetTrajectory(i).back();}
    }
    cout << "Média de 200 partículas ao fim de 10 passos: " << sum/mov1.GetNumber() << endl << endl;

    sum = 0;
    double max2 = 0;
    for(int i = 0; i < mov2.GetNumber(); ++i){
        sum += mov2.GetTrajectory(i).back() - mov2.GetTrajectory(i).front();
        if (mov2.GetTrajectory(i).back() > max2) {max2 = mov2.GetTrajectory(i).back();}
    }
    cout << "Média de 200 partículas ao fim de 100 passos: " << sum/mov2.GetNumber() << endl << endl;

    // ---------------------------
    TCanvas canvas;
    cout << "Histogramas:" << endl << endl;

    TH1D histmov1("posicoes", "posicoes de 200 particulas ao fim de 10 passos;posicao;ocorrencias", 20, 0.0, max1);
    for(int i = 0; i < mov1.GetNumber(); ++i){
        histmov1.Fill(mov1.GetTrajectory(i).back());
    }
    histmov1.SetLineColor(kOrange + 8);
    histmov1.SetFillColor(kOrange + 8);
    histmov1.SetLineWidth(3);
    histmov1.Draw();
    canvas.Update();
    canvas.SaveAs("histmov1.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();

    TH1D histmov2("posicoes", "posicoes de 200 particulas ao fim de 100 passos;posicao;ocorrencias", 20, 0.0, max2);
    for(int i = 1; i <= mov2.GetNumber(); ++i){
        histmov2.Fill(mov2.GetTrajectory(i).back());
    }
    histmov2.SetLineColor(kOrange + 8);
    histmov2.SetFillColor(kOrange + 8);
    histmov2.SetLineWidth(3);
    histmov2.Draw();
    canvas.Update();
    canvas.SaveAs("histmov2.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();
    return 0;
}