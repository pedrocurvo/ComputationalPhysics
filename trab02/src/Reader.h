#ifndef __READER_H__
#define __READER_H__

#include <Eigen/Dense>
#include <string>
#include <vector>
using namespace std;

class Reader{

public:
    // constructor e desctructor:
    Reader(string name, int c): fname(name), cols(c) {;}
    ~Reader() = default;

    // Attributes:

    /*
    Esta função tem como objetivo ler o ficeiro covariance_data.dat e colocar os seus dados num
    vetor de vetores:
    */
    vector<vector<double>> GetData();

protected:

    std::string fname;
    int cols;
};

#endif 