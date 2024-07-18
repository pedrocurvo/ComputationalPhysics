#include <IntegDeriv.h>
#include <iostream>
#include <TRandom.h>
#include <random>
#include <cmath>
#include <iomanip>
using namespace std;

IntegDeriv::IntegDeriv(Functor& f) : F(f) {};

double IntegDeriv::SecondDerivative(double x, double h) {
    return (-F(x-2*h) + 16*F(x-h) - 30*F(x) + 16*F(x+h) - F(x+2*h))/(12*h*h);
}

double IntegDeriv::FirstDerivative(double x, double h){
	return  (F(x - 2*h)+8*F(x + h) - (8*F(x - h) + F(x + 2*h))) / (12*h);
}

// double IntegDeriv::ThirdDerivative(double x, double h){
//     return 0;
// }

double IntegDeriv::FourthDerivative(double x, double h){
    return (F(x + h) + F(x - h) - 2 * F(x) - h*h*IntegDeriv::SecondDerivative(x, h)) * 24 / (2 * pow(h, 4));
}


void IntegDeriv::TrapezoidalRule(double xi, double xf, double& Integral, double& Error, bool control){
    Integral = (xf - xi) / 2 * (F(xi) + F(xf));
        int counter = 1;
    if(control){   
        while(true){
            double h = (xf - xi) / pow(2, counter);
            double sum = 0;
            for(int i = 0; i < pow(2, counter); i++){
                sum += F(xi + i*h) + F(xi + i*h + h);
            }
            sum = sum / 2 * h;
            Integral = sum;
            double erro = 0;
            for(int i = 1; i < pow(2, counter); i++){
                erro += - pow(h, 3) / 12 * IntegDeriv::SecondDerivative(xi + i * h + h/2);
            }
            if(abs(erro) < Error){
                Error = abs(erro);
                break;
            }
            counter++;
        }
    }else{
        while(true){
            double h = (xf - xi) / counter;
            double sum = 0;
            for(int i = 0; i < counter; i++){
                sum += F(xi + i*h) + F(xi + i*h + h);
            }
            sum = sum / 2 * h;
            Integral = sum;
            double erro = 0;
            for(int i = 1; i <= counter; i++){
                erro += - pow(h, 3) / 12  * IntegDeriv::SecondDerivative(xi + i * h + h/2);
            }
            if(abs(erro) < Error){
                Error = abs(erro);
                break;
            }
            counter +=2 ;
        }
    }
}

 void IntegDeriv::simpsonRule(double xi, double xf, double& Integral, double& Error){
    // double h = (xf - xi) / 2;
    // Integral = (F(xi) + 4*F(xi + h) + F(xi + 2*h)) * h / 3;



    int counter = 2;
    while(true){
        double h = (xf - xi) / counter;
        double sum = 0;
        for(int i = 0; i < counter - 2; i+= 2){
            sum += (F(xi +i*h) + 4*F(xi + h+i*h) + F(xi + 2*h + i*h)) ;
        }
        sum = sum * h / 3;
        Integral = sum;
        double erro = 0;
        for(int i = 0; i < counter; i++){
            erro += (xf - xi) * pow(h, 4) / 180 * IntegDeriv::FourthDerivative(xi + i * h + h/2);
        }
        cout << "erro " << erro << endl;
        if(abs(erro) < Error){
            Error = erro;
            break;
        }
        counter += 2;
    }
 }

void IntegDeriv::MonteCarloNormal(double xi, double xf, double& Integral, double& Error){
    int N = 2;
    while(true){
        double sum = 0;
        double erro = 0;
        double square_of_random = 0;
        for(int j = 1; j <= N; ++j){
            double n = gRandom->Rndm();
            n = xi + n * (xf - xi);
            sum += F(n);
            square_of_random += pow(F(n), 2);
        }
        erro = (xf - xi) / sqrt(N) * sqrt(square_of_random / N - pow(sum/N, 2));
        sum *= (xf - xi) / N;
        Integral = sum;
        if(abs(erro) < Error){
            Error = erro;
            break;
        }
        N *= 2;
    }
}

void IntegDeriv::MonteCarloVonNeumann(double xi, double xf, double& Integral, double& Error, int iterations){
    double h = (xf - xi) / abs(2*(int)(xf - xi));
    double max = 0;
    for(int i = 0; i < abs(2*(int)(xf - xi)); i++){
        if(F(xi + h * i) > max){
            max = F(xi + h * i);
        }
    }
    double Area = (xf - xi) * max;
    auto const hes = std::random_device{}();
    gRandom->SetSeed(hes);
    int counter = 0;
    for(int i = 0; i < iterations; i++){
        // double n = gRandom->Rndm();
        // double n_y = gRandom->Rndm();
        double x_random = xi + gRandom->Rndm() * (xf - xi);
        double y_random = max * gRandom->Rndm();
        if(y_random <= F(x_random)){
            counter++;
        };
    }
    Integral = Area * counter / iterations;
    Error = (xf-xi) * max / iterations * sqrt(counter * (1 - (counter/iterations)));
}

