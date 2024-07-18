#include <cmath>
#include <iostream>
#include <string>
#include <Functor.h>
#include <MyFunction.h>
using namespace std;

MyFunction::MyFunction(std::string s, std::function<double(double)> fv) : Functor(s) {
    function = fv;
}

double MyFunction::operator() (double x) {
    return function(x);
}