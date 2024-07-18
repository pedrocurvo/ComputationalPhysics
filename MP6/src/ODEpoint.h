#ifndef __ODEPOINT_H__
#define __ODEPOINT_H__

#include <Xvar.h>
#include <ODEpoint.h>
#include <iostream>
using namespace std;

class ODEpoint : public Xvar {
    
    public:
        ODEpoint() : t(-1) {;}
        ODEpoint(double t_, Xvar a_) : t(t_), Xvar(a_) {;}
        ODEpoint(double t_, const std::vector<double>& v) : t(t_), Xvar(v) {;}
        ODEpoint(double t_, const std::initializer_list<double>& v) : t(t_), Xvar(v) {;}
    
    void SetODEpoint(double t_, Xvar& p);
    void SetODEpoint(double t_, const std::initializer_list<double>& v);
    void SetODEpoint(double t_, std::vector<double> v);
    double& operator[](int); // X[i]
    friend std::ostream& operator<< (std::ostream&, const ODEpoint&);
    
    double& T() { return t;} // accessor to time

    private:
        double t; // time
};

#endif