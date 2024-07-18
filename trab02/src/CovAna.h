#ifndef __COVANA__
#define __COVANA__

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <TGraph2D.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TAxis.h"
#include "TApplication.h"
#include "Reader.h"
#include "TSystem.h"
#include <TStyle.h>
using namespace std;

class CovAna
{
    public:
        // Constructor:
        CovAna() = default;
        CovAna(
                vector <vector <double>>& );

        // Destructor:
        ~CovAna() = default;

        /*
        Retorna a matriz de covariâncias
        */
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> GetCovarianceMatrix();

        /*
        Retorna a média para a variável passada como argumento
        */
        double GetMean(double);

        /*
        Recebe um gráfico e adiciona os pontos ao mesmo
        */ 
        void GetGraph2D(TGraph2D& );
        void DrawHistogram(int index, TH1D& hist, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> R);
    
    private:
        vector <vector <double>> mat;
};

#endif