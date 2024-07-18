#include <Xvar.h>

using namespace std;

Xvar::Xvar(int n){
    vector<double> vec(n);
    x = vec;
};

Xvar::Xvar(const Xvar& original): Xvar(original.x){;}; //perguntar orcinha
Xvar& Xvar::operator=(const Xvar& original){
    this->x = original.x;
    return *this;
}; // assignment operator

Xvar Xvar::operator+(const Xvar& original){ // perguntar orcinha
    Xvar result(x.size());
    for(int i = 0; i < original.x.size(); i++){
        result.x[i] = this->x[i] + original.x[i];
    } //mudar original??
    return result;
}

double& Xvar::operator[](int index){
    return this->x[index];
}

Xvar operator*(double scalar, const Xvar& obj){
    Xvar result(obj.x.size());
    for(int i = 0; i < obj.x.size(); i++){
        result.x[i] = scalar * obj.x[i];
    }
    return result;
}

std::ostream& operator<<(std::ostream&, const Xvar& original){
    cout << "[";
    for(int i = 0; i < original.x.size() - 1; i++){
        cout << original.x[i] << ", ";
    }
    cout << original.x[original.x.size() - 1] << "]"<< endl;
    return cout;
}

std::vector<double>& Xvar::X(){
    return this->x;
}