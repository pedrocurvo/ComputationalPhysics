#include <iostream>
#include <IntegDeriv.h>
#include <cmath>
#include <iomanip>
#include <random>
#include <chrono>
#include <TH1D.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <ctime>
#include <TApplication.h>
#include <TF1.h>
#define DRAW
using namespace std;

int main(){

    //define function to integrate
    auto f = [](double *x, double *par){
        //par[0]; mean
        //par[1]; sigma
        return 1/(sqrt(2*M_PI))/(par[1]) * exp(-0.5*pow((x[0]-par[0]), 2)/(par[1]*par[1]));
    };
    auto F = new TF1("Gaussiana;x;F(x)", f, -5, 5, 2);
    #ifdef DRAW
        TApplication A("A", nullptr, nullptr);
        TCanvas c("canvas", "", 1200, 900);
        F->SetLineWidth(4);
        F->SetLineColor(kRed + 2);
        F->SetFillColor(39);
        F->SetFillStyle(1001);
        //mean is 0
        //sigma is 1
        F->SetParameters(0, 1);
        F->Draw("LF2");
        c.Update();
        c.WaitPrimitive();
    #endif

    //integration limits
    double x1 = -3.0;
    double x2 = 3.0;

    //gaussian params
    double par[2] {0, 1};
    //Number of randoms
    int N = 1e6;
    gRandom->SetSeed(time(NULL));

    //accumulated values of x
    double ftot = 0;
    double ftot2 = 0;
    for(int i = 0; i < N; i++){
        //random number between 0 and 1
        double r = gRandom->Uniform();

        //random number between x1 and x2
        double x = x1 + r * (x2 - x1);

        ftot = ftot + f(&x, par);
        ftot2 = ftot2 + f(&x, par) * f(&x, par);
    }

    //mean value of f
    double fmean = ftot / N;
    double integral = fmean * (x2 - x1);



    cout << "Integral MC " << integral << endl;

    //variance of f = <f*2> - <f>^2
    double var_f = ftot2/N - fmean * fmean;
    double sigma_f_mean = sqrt(var_f/N);

    //Integral error
    double integral_error = (x2 - x1) * sigma_f_mean;
    cout << "Integral error: " << integral_error << endl;
    cout << "Integral: " << integral << " +/- " << integral_error << endl;


    //Generating Random numbers
    TH1D h("h", "Random Numbers", 100, 0, 20);
    // Formas de ter seeds boas
    auto const hes = std::random_device{}();
    auto const seed = std::chrono::system_clock::now().time_since_epoch().count();
    //
    gRandom->SetSeed(hes);
    for(int i = 0; i < N; i++){
        double r = gRandom->Rndm();

        //determine the x corresponding to integral of f(x)
        double x = -log(1 - r);
        h.Fill(x);
    }
    #ifdef DRAW
        gPad->SetLogy();
        h.SetLineColor(kRed + 2);
        h.SetFillColor(39);
        h.SetFillStyle(1001);
        //mean is 0
        //sigma is 1
        h.Draw();
        c.Update();
        c.WaitPrimitive();
    #endif

    return 0;
}