#ifndef __INTERPOLATOR_H__
#define __INTERPOLATOR_H__

#include <iostream>
#include <DataPoints.h>
#include <map>
#include <vector>
#include <initializer_list>
#include <string>
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TSystem.h"
#include <Eigen/Dense>
using namespace std;

class Interpolator : public DataPoints {
    public:
    // constructors, destructor
    Interpolator() = default; //default constructor (nothing to be done?)
    Interpolator(initializer_list<pair<double, double>>);
    Interpolator(int N, double* x, double* y); // build DataPoints from C-arrays of x and y values
    Interpolator(const std::vector< std::pair<double,double> >&);
    Interpolator(const std::vector< double>& x, const std::vector< double>& y);
    ~Interpolator() = default;
    // interpolation methods
    void Init(); // calculations need at init time (a coeffs, map)
    double InterpolateLagrange(double); // Lagrange interpolation
    double InterpolateNewton(double); // Newton interpolation
    double InterpolateSpline3(double); // spline3
    double CosineInterpolate(double mu); // cosine interpolation
    double Lagrange(double *x, double *p);
    double Newton(double *x, double *p);
    double Spline3(double *x, double *p);
    double Cosine(double *x, double *p);


    // // draw points and function
    void Draw(std::string s); // s="lagrange", "newton", "spline3"
    
    
    private:
    std::vector<double> x,y; // data points x,y
    std::map<std::string,TF1*> MI; // key="lagrange", "newton", "spline3"
    // Newton interpolation
    std::vector<double> a; // coefficients
};
 
 #endif