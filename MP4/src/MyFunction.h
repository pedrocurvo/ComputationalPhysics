#ifndef __MYFUNCTION_H__
#define __MYFUNCTION_H__

#include "Functor.h"
#include <functional>
using namespace std;

class MyFunction : public Functor {
    public:
        MyFunction(std::string s, std::function<double(double)>);// : Functor(s) {;}
        ~MyFunction() = default;

        // redefine function
        virtual double operator()(double);
    private:
        std::function<double(double)> function;
};

#endif