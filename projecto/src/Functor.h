#ifndef __FUNCTOR_H__
#define __FUNCTOR_H__

#include <string>
#include "TAxis.h"
#include <iostream>
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraph.h"
using namespace std;

class Functor {
    public:
        Functor(std::string s="Functor") : name(s) {;}
        ~Functor() = default;
        virtual double operator()(double x);
        virtual void Draw(double xi, double xf, int num, std::string xtitle="x", std::string ytitle="y");
    protected:
        static TCanvas *c;
        std::string name;
};

#endif