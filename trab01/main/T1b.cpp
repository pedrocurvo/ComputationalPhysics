#include <iostream>
#include <T1func.h>
#include <numeric>
#include <vector>

//Root includes
#include "TCanvas.h"
#include "TH1D.h"
#include "TApplication.h" 
#include "TGraph.h" 
#include "TSystem.h"
#include "TF1.h"
using namespace std;

int main(){

    // Read data file to arrays
    const int N = 5000;
    double t[N], A[N];
    string fname = "tDFT_sinal.dat"; // file name
    ReadFile(fname, t, A);

    //Get Fourier Coefficients
    vector< pair<double, complex<double>>> vC(N/2);
    GetFourierCoeffs(N, t, A, vC);

    //Reconstruct signal in time, f(t), from detected harmonics
    // using a lambda function and a TF1 
    vector<pair<double, complex<double>>> vCC = important_freq(vC);
    int M {2500};
    auto frec = [&](double* x, double* par){
        complex <double> termo = 0;
        for (int j = 0; j < M; ++j){
            double part_real = cos((2*M_PI)*(vCC[j].first)*x[0]);
            double part_ima = sin((2*M_PI)*(vCC[j].first)*x[0]);
            complex <double> number = part_real + part_ima * 1i;

        termo += (vCC[j].second)*number;
        }
        termo /= M;
        return real(termo);
    };
    
    auto sinal_filtrado = new TF1("Sinal Filtrado",frec,0.0,5.0,0);
    
    //compare data file amplitudes with reconstructed ones 
    cout << " |   time   |  amplitude | signal rec |" << endl;
    for(int i{0}; i < 20; i++){
            cout << setw(12) << t[i] << " " << setw(12) << A[i] << " " << setw(12) <<frec(&t[i], nullptr) << endl;
    }

    ///Drawing
    TApplication App("A",0,0);
	TCanvas canvas("canvas","",1200,800);
    

    TGraph amp_original;
    for (int c = 0; c<5000;c++){
        amp_original.AddPoint(t[c],A[c]);
    }
    amp_original.SetMarkerStyle(25);
    amp_original.SetMarkerSize(1.5);
    amp_original.SetMarkerColor(kBlue+2);
    amp_original.SetLineColor(kBlue+2);
    amp_original.SetLineWidth(3);
    amp_original.Draw("AL");
    amp_original.SetTitle("Amplitude vs Tempo; tempo; amplitude");

    sinal_filtrado->SetLineWidth(4);
    sinal_filtrado->SetLineColor(kRed +2);
    sinal_filtrado->Draw("same");
    canvas.Update();
    canvas.SaveAs("T1b_sinalrec.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();



}







