#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>
#include <vector>
using namespace std;

void ReadFile(string fname, double (&t)[5000], double (&A)[5000]){
    fstream file;
    file.open(fname, fstream::in);
    char line[100];

    file.getline(line, 100);
    file.getline(line, 100);

    int  i{};
    double column1, column2;
    while(file >> column1 >> column2){
        t[i] = column1;
        A[i] = column2;
        ++i; 
    };
    file.close();
}

void ReadFile(string fname, double (&t)[2000], double (&A)[2000]){
    fstream file;
    file.open(fname, fstream::in);
    char line[100];

    file.getline(line, 100);
    file.getline(line, 100);

    int  i{};
    double column1, column2;
    while(file >> column1 >> column2){
        t[i] = column1;
        A[i] = column2;
        ++i; 
    };
    file.close();
}


vector <double> media_deslizante(int N , double *A, int band = 3){
    band = (int)(band / 2);
    vector <double> media_D(N);
    for(int i {0}; i < N; i++){
        double sum = 0;
        if(i < band || i > N - band){
            media_D[i] = A[i];
        }else{
            for(int j{i - band}; j <= i + band; j++){
                sum += A[j];
            }
            sum /= (1 + 2 * band);
            media_D[i] = sum;
        }
    }
    return media_D;
}

double GetTimeSampling(string fname){
    fstream file;
    file.open(fname, fstream::in);
    double column1, column2, column3, column4;
    file >> column1 >> column2 >> column3 >> column4;
    return column4;
}


double autocorrelacao(int k, double T, double *A){
    double media {}, nom{}, denom{};
    for (int i{}; i < T; i++){
        media += A[i] / T;
    }
    for(int t {k + 1}; t <= T; t++ ){
        nom += (A[t] - media) * (A[t-k] - media);
    }
    for(int t{1}; t < T; t++){
        denom += pow((A[t] - media), 2);
    }
    return (T / (T - k)) * nom / denom;
};

vector <pair<double, double>> GetAutoCorrelation(const int N, double *t, double *A){
    vector <pair<double, double>> vAuto;
    for(int i{}; i < N; i++){
        pair<double, double> ts_ac;
        ts_ac.first = i * 0.001;
        ts_ac.second = autocorrelacao(i, N, A);
        vAuto.push_back(ts_ac);
    }
    return vAuto;
}

void GetFourierCoeffs(int N, double (&t)[5000], double (&A)[5000], vector<pair<double,complex<double>>> &vC){
    double part_real;
    double part_ima;
    double total_time = 5.0;
    double sampling_time = 0.001;
    double fs {1/sampling_time};

    for (int i = 0; i < (N/2); ++i){
        double fk {(i*fs)/N};
        complex <double> fourier_coeff = 0;
        for (int j = 0; j < N; ++j){
            part_real = cos(((-2*M_PI)/N)*j*i);
            part_ima = sin(((-2*M_PI)/N)*j*i);
            
            complex <double> number = part_real + part_ima * 1i;
            fourier_coeff += A[j]*(number);
        }
        vC[i].first = fk;
        vC[i].second = fourier_coeff;
    }
}

void GetFourierCoeffs(int N, double (&t)[2000], double (&A)[2000], vector<pair<double,complex<double>>> &vC){
    double part_real;
    double part_ima;
    double total_time = 5.0;
    double sampling_time = 0.001;
    double fs {1/sampling_time};

    for (int i = 0; i < (N/2); ++i){
        double fk {(i*fs)/N};
        complex <double> fourier_coeff = 0;
        for (int j = 0; j < N; ++j){
            part_real = cos(((-2*M_PI)/N)*j*i);
            part_ima = sin(((-2*M_PI)/N)*j*i);
            
            complex <double> number = part_real + part_ima * 1i;
            fourier_coeff += A[j]*(number);
        }
        vC[i].first = fk;
        vC[i].second = fourier_coeff;
    }
}

vector <pair<double,double>> periodograma(double(&A)[5000]){

	vector <pair<double,double>> vP;
	for (int k = 0; k<2500;k++){
		double soma_real = 0;
		double soma_complexa = 0;
		for (int n = 0; n< 5000; n++){
			soma_real += A[n]*cos(-2*(M_PI/5000)*n*k);
			soma_complexa += A[n]*sin(-2*(M_PI/5000)*n*k);
		}
		double fk = (k*1000)/5000;
		//fs = 1000 e N=5000
		double p = (soma_real*soma_real + soma_complexa*soma_complexa)/5000;
		vP.push_back(make_pair(fk,p));
	}
    return vP;
};

vector <pair<double,double>> periodograma(double(&A)[2000]){

    vector <pair<double,double>> vP;
    for (int k = 0; k<1000;k++){
        double soma_real = 0;
        double soma_complexa = 0;
        for (int n = 0; n< 2000; n++){
            soma_real += A[n]*cos(-2*(M_PI/2000)*n*k);
            soma_complexa += A[n]*sin(-2*(M_PI/2000)*n*k);
        }
        double fk = (k*400)/2000;
        double p = (soma_real*soma_real + soma_complexa*soma_complexa)/2000;
        vP.push_back(make_pair(fk,p));
    }
    return vP;
};

bool sort_complex(complex<double> a, complex<double> b){
    //sqrt( pow())
    return norm(a) > norm(b);
}
bool sort_complex2(pair<double, complex<double>> a, pair<double, complex<double>> b){
    //sqrt( pow())
    return norm(a.second) > norm(b.second);
}

vector<complex<double>> important_coeff (vector<pair<double, complex<double>>> &vC){
        vector< complex<double> > res;
        for(int i {0}; i < vC.size(); i++){ 
            res.push_back(vC[i].second);
        };
        sort(res.begin(), res.end(), sort_complex);
        return res;
};

vector<pair<double, complex<double>>> important_freq (vector<pair<double, complex<double>>> &vC){
    vector<pair<double, complex<double>>> res;
    for(int i {0}; i < vC.size(); i++){ 
            res.push_back(vC[i]);
        };
    sort(res.begin(), res.end(), sort_complex2);
    return res;
};




















