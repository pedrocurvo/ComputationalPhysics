#include <Xvar.h>
#include <ODEpoint.h>
#include <iostream>
using namespace std;

void ODEpoint::SetODEpoint(double t_, Xvar& p) {
    t = t_;
    this->x = p.X();
    // Xvar::operator=(p);
}

void ODEpoint::SetODEpoint(double t_, const std::initializer_list<double>& v){
    t = t_;
    this->x = v;
    // Xvar::operator=(v);
}

void ODEpoint::SetODEpoint(double t_, std::vector<double> v){
    t = t_;
    Xvar::operator=(v);
}

std::ostream& operator<<(std::ostream&, const ODEpoint& original){
    cout << "[";
    cout << original.t << ": ";
    for(int i = 0; i < original.x.size() - 1; i++){
        cout << original.x[i] << ", ";
    }
    cout << original.x[original.x.size() - 1] << "]"<< endl;
    return cout;
}

double& ODEpoint::operator[](int index){
    if(index == 0){
        return t;
    }else{
        return this->X()[index];
    };
}

