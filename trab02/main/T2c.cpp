#include <string>
#include <vector>
#include <Eigen/Dense>
#include <FCmatrixAlgo.h>
#include <CovAna.h>
#include <iostream>
#include <TRandom.h>
#include <TStyle.h>
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TApplication.h"
#include "Reader.h"
#include "TSystem.h"
#include "RootTools.h"
#include "Math/GenVector/Cartesian3D.h"
using namespace std;

int main(){
    string filename = "covariance_data.dat";
    Reader RR(filename, 3); //Read file with 3 columns per line 
    vector<vector<double>> data = RR.GetData();

    //----------------------------------------------------------------
    //Produce Covariance Matrix
    //----------------------------------------------------------------
    CovAna COV(data);
    auto V = COV.GetCovarianceMatrix();
    // print covariance matrix
    cout << "Covariance Matrix: " << endl << V << endl << endl;
    //----------------------------------------------------------------
    //print variable means (medias)
    cout << "<x1>=" << COV.GetMean(1) << endl;
    cout << "<x2>=" << COV.GetMean(2) << endl;
    cout << "<x3>=" << COV.GetMean(3) << endl;
    cout << endl;

    //----------------------------------------------------------------
    //Eigen values and Eigenvectors analysis
    //----------------------------------------------------------------
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Q(V.rows(), V.cols());
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> R(V.rows(), V.cols());
    FCmatrixAlgo::QRfactorization(V, Q, R);
    // factorize V matrix in QR and print results
    cout << "V: " << endl << V << endl << endl;
    cout << "Q: " << endl << Q << endl << endl;
    cout << "R: " << endl << R << endl << endl;
    //Show that is correct QR
    cout << "Q*R(should be equal to V): " << endl << Q*R << endl << endl;
    //----------------------------------------------------------------
    //Show that Q is orthogonal
    //----------------------------------------------------------------
    cout << "Confirming that is orthogonal: " << endl;
    cout << "Q*QT (should be equal to identity (might not have zeros due to CPU precision)): " << endl << Q.transpose() * Q << endl << endl;
    //----------------------------------------------------------------
    //Get Eigenvalues and Eigenvectors
    //----------------------------------------------------------------
    double eps = 1;
    FCmatrixAlgo::IterativeQR(V, Q, R, eps);
    //Eigenvalues
    cout << "------------Eigenvalues--------------------" << endl;
    for(int i = 0; i < Q.rows(); i++){
        cout << "Eigenvalue " << i + 1 << ": " << Q(i,i) << endl;
    }
    cout << endl;

    cout << "------------Eigenvectors------------------" << endl;
    //Eigenvectors
    for(int row = 0; row < R.cols(); row++){
        cout << "Eigenvector " << row + 1 << ": " << endl << R.col(row) << endl;
    }

    TApplication A("k", nullptr, nullptr);
    
    TCanvas *c = new TCanvas("c","Data Points",0,0,1200,900);
    TGraph2D *G = new TGraph2D();
    COV.GetGraph2D(*G);
    G->Draw("pcol");
    CanvasUpdater(c, "T2c_1.pdf");
    //----------------------------------------------------------------
    // Draw vectors
    //----------------------------------------------------------------
    //One of the vectors is irrelevant due to normalization by it's eigenvalue
    for(int i = 0; i < R.cols(); i++){
        double x = R(0,i);
        double y = R(1,i);
        double z = R(2,i);
        TDrawVector(x * Q(i, i), y * Q(i, i), z * Q(i, i));
    }
    G->SetTitle("Data Points with Eigenvectors");
    CanvasUpdater(c, "T2c_3.pdf");
    //----------------------------------------------------------------
    TH1D hist1("hist1", "Distancias 1", 100, 0, 100);
    COV.DrawHistogram(0, hist1, R);
    CanvasUpdater(c, "T2c_21.pdf");
    //----------------------------------------------------------------
    TH1D hist2("hist2", "Distancias 2", 100, 0, 100);
    COV.DrawHistogram(1, hist2, R);
    CanvasUpdater(c, "T2c_22.pdf");
    //----------------------------------------------------------------
    TH1D hist3("hist3", "Distancias 3", 100, 0, 100);
    COV.DrawHistogram(2, hist3, R);
    CanvasUpdater(c, "T2c_23.pdf");

    //--------------------------Important Directions-------------------------------------
    TGraph *G2 = new TGraph();
    vector <vector <double>> pontos2D;
    /*
    Vamos agora escrever as coordenadas dos nossos dados na base dos vetores próprios da matriz de
    covariâncias. Deste modo, podemos desenhar gráficos que reflitam evidentemente a variação
    dos dados nestas direções:
    */
    for(int i = 0; i < data.size(); ++i)
    {
        Eigen::Matrix<double,3,1> dados = {data[i][0], data[i][1], data[i][2]};
        Eigen::Matrix<double,3,1> res = R.transpose() * dados;
        vector <double> vetor = {res[0], res[1]};
        pontos2D.push_back(vetor);
    }

    for (auto& values: pontos2D)
    {
        G2->AddPoint(values[0], values[1]);
    }
    G2->SetTitle("Graph Data in 2D; X (direction in one E.V.); Y (direction in other E.V.)");
    G2->SetMarkerStyle(4);
    G2->SetMarkerColor(kBlue);
    G2->SetMarkerSize(1.);
    G2->Draw("AP");
    CanvasUpdater(c, "T2c_41.pdf");
    gSystem->ProcessEvents();
    //-----------------------------------------------------------------------------------
    //-------------------------- Directions that doesnt matter-------------------------------------
    TGraph *G3 = new TGraph();
    vector <vector <double>> pontos2D1;
    /*
    Vamos agora escrever as coordenadas dos nossos dados na base dos vetores próprios da matriz de
    covariâncias. Deste modo, podemos desenhar gráficos que reflitam evidentemente a variação
    dos dados nestas direções:
    */
    for(int i = 0; i < data.size(); ++i)
    {
        Eigen::Matrix<double,3,1> dados = {data[i][0], data[i][1], data[i][2]};
        Eigen::Matrix<double,3,1> res = R.transpose() * dados;
        vector <double> vetor = {res[1], res[2]};
        pontos2D1.push_back(vetor);
    }

    for (auto& values: pontos2D1)
    {
        G3->AddPoint(values[0], values[1]);
    }
    G3->SetTitle("Graph Data in 2D (irrelevant direction); X (direction in one E.V.); Y (direction in other E.V.)");
    G3->SetMarkerStyle(4);
    G3->SetMarkerColor(kRed);
    G3->SetMarkerSize(1.);
    G3->Draw("AP");
    CanvasUpdater(c, "T2c_42.pdf");
    gSystem->ProcessEvents();
    //-----------------------------------------------------------------------------------
    TGraph *G4 = new TGraph();
    vector <vector <double>> pontos2D2;
    /*
    Vamos agora escrever as coordenadas dos nossos dados na base dos vetores próprios da matriz de
    covariâncias. Deste modo, podemos desenhar gráficos que reflitam evidentemente a variação
    dos dados nestas direções:
    */
    for(int i = 0; i < data.size(); ++i)
    {
        Eigen::Matrix<double,3,1> dados = {data[i][0], data[i][1], data[i][2]};
        Eigen::Matrix<double,3,1> res = R.transpose() * dados;
        vector <double> vetor = {res[0], res[2]};
        pontos2D2.push_back(vetor);
    }

    for (auto& values: pontos2D2)
    {
        G4->AddPoint(values[0], values[1]);
    }
    G4->SetTitle("Graph Data in 2D (irrelevant direction); X (direction in one E.V.); Y (direction in other E.V.)");
    G4->SetMarkerStyle(4);
    G4->SetMarkerColor(kGreen +3);
    G4->SetMarkerSize(1.);
    G4->Draw("AP");
    CanvasUpdater(c, "T2c_43.pdf");
    gSystem->ProcessEvents();

    A.Run();
    return 0;
}