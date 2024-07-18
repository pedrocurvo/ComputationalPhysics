#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "EqSolver.h"
#include "FCmatrixAlgo.h"
using namespace std;

int main() {
    // Creation of matrix and vector of constants
    int dim = 4;
    matriz_din A(dim, dim);
    vetor_coluna b(dim);
    double V0 = 20;
    double R1 = 8.73366;
    double R2 = 13.0173;
    double R3 = 10.239;
    double R4 = 3.58372;
    double R5 = 1.81721;
    double R6 = 9.78139;
    double R7 = 9.26458;
    double R8 = 43.5632;
    // Fill matrix and vector
    A << R1 + R2, -R1, -R2, 0,
         -R1 , R1 + R3 + R6 + R4, -R4, -R6,
         -R2, -R4, R2 + R4 + R7 + R5, -R7,
        0, -R6, -R7, R8 + R7;
    b << V0, 0, 0, 0;
    // Object Solver 
    EqSolver S(A, b);
    //------------------------------
    // Get solution by inverse
    //------------------------------
    matriz_din A_inv(dim, dim);
    FCmatrixAlgo::Invert(A, A_inv);
    auto x1 = A_inv * b;
    cout << "Solution by inverse: " << endl;
    cout << x1 << endl << endl;
    cout << A_inv  << endl;
    cout << A_inv * A << endl;

    return 0;
}