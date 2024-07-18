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
    ReadFile("tDFT_sinal.dat", t, A);
    vector<double> media_des = media_deslizante(N, A, 21);
    vector<pair<double, double>> vA = GetAutoCorrelation(N, t, A);

    cout << "20 1st values of time and auto_correlation" << endl;
    for(int i = 0; i < 20; i++){ 
        cout << "Time: " << setw(7) << vA[i].first << "   Auto-correlation : "  << setw(10) << vA[i].second << endl;
    }
	cout << endl;

    //Get Fourier Coefficients
    vector<pair<double, complex<double>>> vC(N/2); 
    GetFourierCoeffs(N, t, A, vC);
  
	//Print 6 most important values of Fourier Coeficientes (the ones with highest absolute value)
    vector < complex<double> > six_cf = important_coeff(vC);
    for_each(six_cf.begin(), six_cf.begin() + 6, [](complex<double> x){cout << "F.Coef : " << x << endl;});


    vector <pair<double,double>> vP = periodograma(A);

    TApplication App("A",0,0);
	TCanvas canvas("canvas","",1200,800);

	//1st graph 
	TGraph autocorrelacao_vs_tempo;
	for(int i = 0 ;i < 200; i++){
		autocorrelacao_vs_tempo.AddPoint(t[i], vA[i].second);
	}
	autocorrelacao_vs_tempo.SetMarkerStyle(25);
    autocorrelacao_vs_tempo.SetMarkerSize(1.5);
    autocorrelacao_vs_tempo.SetMarkerColor(kGreen+2);
    autocorrelacao_vs_tempo.SetLineColor(kGreen+2);
    autocorrelacao_vs_tempo.SetLineWidth(3);
    autocorrelacao_vs_tempo.Draw("AL");
    autocorrelacao_vs_tempo.SetTitle("Autocorrelacao vs Tempo; tempo; autocorrelacao");
    canvas.Update();
    canvas.SaveAs("T1a_autocorrelacao.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();

    // 2nd graph
	TGraph sinal_vs_tempo;
	for(int i = 0; i < 200; ++i){
	    sinal_vs_tempo.AddPoint(t[i], A[i]);
	}
	sinal_vs_tempo.SetMarkerStyle(25);
	sinal_vs_tempo.SetMarkerSize(1.5);
	sinal_vs_tempo.SetMarkerColor(kBlue+2);
	sinal_vs_tempo.SetLineColor(kBlue+2);
	sinal_vs_tempo.SetLineWidth(3);
	sinal_vs_tempo.Draw("APL");
	sinal_vs_tempo.SetTitle("Tempo vs Amplitude; tempo; amplitude");

	TGraph media_deslizante_graph;
	for (int a = 0; a < 200 ; a++){
		media_deslizante_graph.AddPoint( t[a] , media_des[a]);
	}
	media_deslizante_graph.SetMarkerStyle(25);
	media_deslizante_graph.SetMarkerSize(1.5);
	media_deslizante_graph.SetMarkerColor(kRed+2);
	media_deslizante_graph.SetLineColor(kRed+2);
	media_deslizante_graph.SetLineWidth(3);
	media_deslizante_graph.Draw("same"); // axis drawn
	media_deslizante_graph.SetTitle("Tempo vs Media Deslizante; tempo; media deslizante");
	canvas.Update();
	canvas.SaveAs("T1_a_sinal.png");
	canvas.WaitPrimitive();
	gSystem->ProcessEvents();

	// periodogramas vs frequência
	TGraph periodograma_vs_freq;

	for(int b = 0; b < vC.size(); ++b){
	    periodograma_vs_freq.AddPoint(vP[b].first, vP[b].second);
	}
    //canvas.SetLogy();
	periodograma_vs_freq.SetMarkerStyle(25);
	periodograma_vs_freq.SetMarkerSize(1.5);
	periodograma_vs_freq.SetMarkerColor(kBlue+2);
	periodograma_vs_freq.SetLineColor(kBlue+2);
	periodograma_vs_freq.SetLineWidth(3);
	periodograma_vs_freq.Draw("APL"); // axis drawn
	periodograma_vs_freq.SetTitle("Periodograma vs Freq; frequencia (fk); Densidade Espetral");
	canvas.Update();
	canvas.SaveAs("T1a_periodograma.png");
	canvas.WaitPrimitive();
	gSystem->ProcessEvents();

    return 0;
}