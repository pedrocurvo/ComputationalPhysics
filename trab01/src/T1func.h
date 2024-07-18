#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <complex>
#include <numeric>
#include <vector>
using namespace std;

void ReadFile(string fname, double (&t)[5000], double (&A)[5000]);
void ReadFile(string fname, double (&t)[2000], double (&A)[2000]);
void GetFourierCoeffs(int N, double (&t)[5000], double (&A)[5000], vector<pair<double,complex<double>>> &vC);
void GetFourierCoeffs(int N, double (&t)[2000], double (&A)[2000], vector<pair<double,complex<double>>> &vC);

double GetTimeSampling(string fname);
vector <double> media_deslizante(int N , double *A, int band = 3);
vector <pair<double, double>> GetAutoCorrelation(const int N, double *t, double *A);
vector <pair<double,double>> periodograma(double(&A)[5000]);
vector <pair<double,double>> periodograma(double(&A)[2000]);
vector<complex<double>> important_coeff (vector<pair<double, complex<double>>> &vC);
bool sortbyfirst(const pair<double, complex<double>> &a, const pair<double, complex<double>> &b);
vector<pair<double, complex<double>>> important_freq (vector<pair<double, complex<double>>> &vC);