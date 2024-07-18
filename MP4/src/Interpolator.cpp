#include <Interpolator.h>
using namespace std;

Interpolator::Interpolator(initializer_list< pair<double, double>> L): DataPoints(L){
    Init();
}

Interpolator::Interpolator(int N, double* x, double* y): DataPoints(N, x, y) {
    Init();
}

Interpolator::Interpolator(const std::vector< std::pair<double,double> >& L) : DataPoints(L) {
    Init();
}

Interpolator::Interpolator(const std::vector< double>& x, const std::vector< double>& y) : DataPoints(x, y) {
    Init();
}



void Interpolator::Init(){
    for(const auto&e : P){
        x.push_back(e.first);
        y.push_back(e.second);
    }
}


double Interpolator::InterpolateLagrange(double t){
    double resultado = 0;
    for(int i = 0; i < x.size(); i++){
        double aux_res = 1;
        for(int j = 0; j < y.size(); j++){
            if(i != j){
                aux_res *= (t - x[j])/(x[i] - x[j]);
            }
        }
        resultado += aux_res * y[i];
    }
    return resultado;
}


double Interpolator::InterpolateNewton(double z){
    a = y;
    int n = a.size();
    for(int k = 1; k < n + 1; k++){
        for(int i = k; i < n + 1; i++){
            a[i] = (a[i] - a[k - 1]) / (x[i] - x[k - 1]);
        }
    }
    double P = a[n - 1];
    for(int k = 1; k < n+1; k++){
        P = a[n-k] + (z - x[n-k]) * P;
    }
    return P;
}


double Interpolator::InterpolateSpline3(double z){
    int n = x.size() - 2;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A(n, n + 1);
    Eigen::Matrix<double, Eigen::Dynamic, 1> Res(n, 1);
    A(0, 0) = 2 * (x[0] - x[2]);
    A(0, 1) = (x[1] - x[2]);
    for(int row = 1; row < A.rows(); row++){
        A(row, row - 1) = (x[row] - x[row + 1]);
        A(row, row) = 2*(x[row] - x[row + 2]);
        A(row, row + 1) = (x[row + 1] - x[row + 2]);
    }
    A.conservativeResize(n, n);
    for(int i = 1; i <= n ; i++){
        Res(i - 1, 0) = 6 * ( (y[i - 1] - y[i]) / (x[i - 1] - x[i]) - (y[i] - y[i + 1])/(x[i] - x[i + 1]) );
    }
    auto K = A.inverse() * Res;
    Eigen::Matrix<double, Eigen::Dynamic, 1> K_final(x.size(), 1);

    K_final(0) = 0;
    for(int i = 0; i < K.rows() ; i++){
        K_final(i + 1) = K(i);
    }
    K_final(K_final.size() - 1) =  0;


    int index = 0;
    for(int i = 0; i < x.size() - 1; i++){
        if(z >=x[i] && z<= x[i+1]){
            index = i;
            break;
        }
    }
    double first = pow(z - x[index + 1], 3) / (x[index] - x[index + 1]) - (z - x[index + 1]) * (x[index] - x[index + 1]);
    double second = pow(z - x[index], 3) / (x[index] - x[index+1]) - (z - x[index]) * (x[index] - x[index + 1]);
    double third = (y[index] * (z - x[index + 1]) - y[index + 1] * (z -x[index])) / (x[index] - x[index+1]);
    double res = K_final[index] * first / 6 - K_final[index + 1] / 6  * second + third;
    return res;
}

void Interpolator::Draw(string var)
{
    TF1 *lagrange = new TF1("lagrange", this, &Interpolator::Lagrange, x[0], x[x.size()-1], 1);
    TF1 *newton = new TF1("newton", this, &Interpolator::Newton, x[0], x[x.size()-1], 2);
    TF1 *spline3 = new TF1("spline3", this, &Interpolator::Spline3, x[0], x[x.size()-1], 2);
    MI.insert({"Lagrange", lagrange});
    MI.insert({"Newton", newton});
    MI.insert({"Spline3", spline3});
    MI["Lagrange"]->SetLineColor(kGreen +3);
    MI["Lagrange"]->SetLineWidth(4);
    MI["Newton"]->SetLineColor(kViolet -7);
    MI["Newton"]->SetLineWidth(4);
    MI["Spline3"]->SetLineColor(kBlue);
    MI["Spline3"]->SetLineWidth(4);

    TGraph g;
    for(int i = 0; i < x.size(); i++)
    {
        g.AddPoint(x[i], y[i]);
    }

    
    g.SetMarkerStyle(25);
    g.SetMarkerSize(2.);
    g.SetMarkerColor(kBlue+2);
    g.SetTitle(var.c_str());
    g.GetXaxis()->SetTitle("x");
    g.GetYaxis()->SetTitle("f(x)");
    g.GetXaxis()->SetTitleSize(0.04);
    g.GetYaxis()->SetTitleSize(0.04);
    g.GetXaxis()->SetTitleOffset(0.6);
    g.GetYaxis()->SetTitleOffset(0.6);
    g.Scale(1/g.GetMaximum());
    TApplication A("Data Points", nullptr, nullptr);
    if(!c) c = new TCanvas("c", "Data Points", 0, 0, 1200, 900);
    g.Draw("AP");
    

    if(var != "Todos"){
        MI[var]->Draw("same");
    }
    else{
        for(auto&e : MI){
            e.second->Draw("same");
        }
    }
    c->Update();
    c->WaitPrimitive();
    gSystem->ProcessEvents();
    A.Run();
}

double Interpolator::Lagrange(double *z, double *p) {
      return this->InterpolateLagrange(*z);
   }
   
double Interpolator::Newton(double *z, double *p) {
      return this->InterpolateNewton(*z);
   }

double Interpolator::Spline3(double *z, double *p) {
      return this->InterpolateSpline3(*z);
   }








