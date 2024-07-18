#include "CovAna.h"
#include "FCmatrixAlgo.h"
#include <iostream>
using namespace std;

CovAna::CovAna(vector<vector<double>>& matrix)
{
    mat = matrix;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> CovAna::GetCovarianceMatrix()
{
    
    int n = mat[0].size();

    // create covariance vector:
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covariante(n,n);
    // create mean vector:
    vector<double> medias(n);
    // Fill mean vector:
    for (int i = 0; i < n; i++)
    {
        medias[i] = CovAna::GetMean(i + 1);
    }
    // Fill covariance matrix:
    for(int row = 0; row < mat[0].size(); ++row)
    {
        for(int col = 0; col < mat[0].size(); ++col)
        {
            double res = 0;
            for(int i = 0; i < mat.size(); ++i)
            {
                res += (mat[i][row] - medias[row])*(mat[i][col] - medias[col]);
            }
            res /= (mat.size() - 1);
            covariante(row,col) = res;
        }
    }
    return covariante;
}

double CovAna::GetMean(double coluna)
{
    double media = 0;
    for(int row = 0; row < mat.size(); ++row)
    {
        media += (mat[row][coluna] / mat.size());
    }
    return media;
}

void CovAna::GetGraph2D(TGraph2D& graph)
{
    // Fill graph:
    for (auto& values: mat)
    {
        graph.AddPoint(values[0], values[1], values[2]);
    }
    graph.SetTitle("Data Points; X; Y; Z");
    graph.SetMarkerStyle(4);
    graph.SetMarkerSize(1.);
    gStyle->SetPalette(kDeepSea);
}
void CovAna::DrawHistogram(int index, TH1D& hist, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> R){
    for(int i = 0; i < 5000; ++i){
        vetor_coluna vetor_distancia(3);
        vetor_distancia(0,0) = mat[i][0] - GetMean(1);
        vetor_distancia(1,0) = mat[i][1] - GetMean(2);
        vetor_distancia(2,0) = mat[i][2] - GetMean(3);
        double ponto = vetor_distancia.transpose() * R.col(index);
        hist.Fill(ponto);
    }
    hist.GetXaxis()->SetTitle("Distancia");
    hist.GetYaxis()->SetTitle("No Pontos");
    hist.SetFillColor(kViolet - 4);
    hist.SetLineColor(kViolet - 6);
    hist.SetLineWidth(2);
    hist.Draw();
}