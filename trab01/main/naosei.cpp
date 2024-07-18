#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iostream>
#include <cmath>

#include "TCanvas.h"
#include "TH1D.h"
#include "TApplication.h" 
#include "TGraph.h" 
#include "TSystem.h"
#include "TF1.h"

//periodograma devolve pair de doubles 

vector <pair<double,double>> periodograma(double(&A)[5000]){
	vector <pair<double,double>> vP;
	for (int k = 0; k<2500;k++){
		double soma_real = 0;
		double soma_complexa = 0;
		for (int n = 0; n< 5000; n++){
			soma_real += A[k]*cos(-(M_PI/1000)*n*k);
			soma_complexa += amp[k]*sin(-(M_PI/1000)*n*k);
		}
		double fk = (k*1000)/5000;
		//fs = 1000 e N=5000
		double p = (soma_real*soma_real + soma_complexa*soma_complexa)/5000;
		vP.push_back(make_pair(fk,p));
	}
}

vector <double> media_deslizante(, int largura){
    vector <double> novo(ampv.begin(),ampv.begin()+largura-1);
    //criar com os primeiros pontos
    int comp = (int) ampv.size();
    double d = 2*largura +1;
    for (int i = largura; i<comp-largura ; i++){
        double sum = accumulate(ampv.begin()+i-largura,ampv.begin()+i+largura,0.0);
        novo.push_back(sum/d);
    }
    for (int i = comp-largura ; i<comp;i++){
        novo.push_back(ampv[i]);
    }
    cout << novo.size() << endl;
    return novo;
}























