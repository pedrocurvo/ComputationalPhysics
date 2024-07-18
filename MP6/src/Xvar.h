#ifndef __XVAR_H__
#define __XVAR_H__

#include <iostream>
#include <vector>
#include <vector>
using namespace std;

class Xvar {
    public:
        Xvar() = default;
        Xvar(int); // number of dependent variables
        Xvar(std::vector<double> v): x(v) {;}; // passing vector
        // using initializer list to build object: X({1,2})
        Xvar(const std::initializer_list<double>& v): x(v) {;};
        ~Xvar() = default;

        Xvar(const Xvar&); // copy constructor
        Xvar& operator=(const Xvar&); // assignment operator
        Xvar operator+(const Xvar&); // operator+
        double& operator[](int); // X[i]

        friend Xvar operator*(double, const Xvar&); // scalar*X
        friend std::ostream& operator<< (std::ostream&, const Xvar&);

        std::vector<double>& X(); // accessor to x
    protected:
        std::vector<double> x;
};

#endif