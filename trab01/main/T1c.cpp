#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <T1func.h>
using namespace std;
// root includes
#include "TCanvas.h"
#include "TH1D.h"
#include "TApplication.h" 
#include "TGraph.h" 
#include "TSystem.h"
#include "TF1.h"

int main(){

    TApplication App2("B",0,0);
    TCanvas canvas("canvas","",1200,800);

    const int N = 2000;
    double t[N], A[N];
    ReadFile("tDFT_sinalruido.dat", t, A);

    //calcular a media deslizante para varios valores de k
    TGraph sum_K;
    for (int d = 0; d<50;d++){
        vector<double> media_des = media_deslizante(N, A, d);
        double soma = 0;
        for (int e = 0; e<2000; ++e) {
            soma += abs(A[e]-media_des[e]);

        }
        sum_K.AddPoint(d,soma);
    }
    sum_K.SetMarkerStyle(1);
    sum_K.SetMarkerSize(1.5);
    sum_K.SetMarkerColor(kBlue+2);
    sum_K.SetLineColor(kBlue+2);
    sum_K.SetLineWidth(3);
    sum_K.Draw("APC");
    sum_K.SetTitle("variacao da md com k; k; SAD");
    canvas.Update();
    canvas.SaveAs("k_SAD.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();


    //fazer a soma da diferença dos valores do grafico original e o novo grafico ponto a ponto
    //fazer um grafico com a soma em funcao do k- ver para que k e q a curva da soma flattens 

    vector <double> media_des2 = media_deslizante(N, A, 10);
    //usamos 10 porque e o valor a partir do qual a diferença entre a media deslizante e o valor original começa a estabilizar- metodo SAD
    

    //criar um array com os novos valores da amplitude obtidos a partir da
    //media deslizante
    double novo_A[2000];
    for (int b = 0; b<2000 ; b++){
        novo_A[b] = media_des2[b];
    }

    //get fourier coefs
    vector<pair<double, complex<double>>> vC(N/2);
    GetFourierCoeffs(N,t,novo_A,vC);

    //criar vetor do periodograma
    vector <pair<double,double>> vP = periodograma(novo_A);

    //concluir qual e o numero de picos a partir do periodograma
    
    

    TGraph sinalruido_vs_tempo;
    for(int i = 0; i < 2000; ++i){
        sinalruido_vs_tempo.AddPoint(t[i], A[i]);
    }

    sinalruido_vs_tempo.SetMarkerStyle(25);
    sinalruido_vs_tempo.SetMarkerSize(1.5);
    sinalruido_vs_tempo.SetMarkerColor(kBlue+2);
    sinalruido_vs_tempo.SetLineColor(kBlue+2);
    sinalruido_vs_tempo.SetLineWidth(3);
    sinalruido_vs_tempo.Draw("APL");
    sinalruido_vs_tempo.SetTitle("Tempo vs Amplitude; tempo; amplitude");
    canvas.Update();
    canvas.SaveAs("tempo_vs_amplitude.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();


    TGraph media_deslizante_graph2;
    for (int a = 0; a < 2000 ; a++){
        media_deslizante_graph2.AddPoint( t[a] , media_des2[a]);
    }
    media_deslizante_graph2.SetMarkerStyle(25);
    media_deslizante_graph2.SetMarkerSize(1.5);
    media_deslizante_graph2.SetMarkerColor(kRed+2);
    media_deslizante_graph2.SetLineColor(kRed+2);
    media_deslizante_graph2.SetLineWidth(3);
    media_deslizante_graph2.Draw("same"); // axis drawn
    media_deslizante_graph2.SetTitle("Tempo vs Media Deslizante; tempo; media deslizante");
    canvas.Update();
    canvas.SaveAs("T1_a_sinal.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();


    //fazer o grafico do periodograma e ver o numero de picos
    TGraph periodograma_freq;
    for (int c = 0; c<2000;c++){
        periodograma_freq.AddPoint(vP[c].first ,vP[c].second);
    }
    periodograma_freq.SetMarkerStyle(25);
    periodograma_freq.SetMarkerSize(1.5);
    periodograma_freq.SetMarkerColor(kBlue+2);
    periodograma_freq.SetLineColor(kBlue+2);
    periodograma_freq.SetLineWidth(3);
    periodograma_freq.Draw("APL"); // axis drawn
    periodograma_freq.SetTitle("Periodograma vs Freq; frequencia (Fk); Densidade Espetral");
    canvas.Update();
    canvas.SaveAs("T1c_periodograma.png");
    canvas.WaitPrimitive();
    gSystem->ProcessEvents();


}









