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

    //----------------------------------------------------------------
    // Get solution by Gauss Elimination
    //----------------------------------------------------------------
    auto x1 = S.GaussSolver(true);
    cout << "Solution by Gauss Elimination: " << endl;
    cout << x1 << endl << endl;
    //----------------------------------------------------------------
    // Get solution by Gauss Elimination Pivoting
    //----------------------------------------------------------------
    auto x2 = S.GaussSolver(false);
    cout << "Solution by Gauss Elimination Pivoting: " << endl;
    cout << x2 << endl << endl;
    //----------------------------------------------------------------
    // Get solution by LU Decomposition
    //----------------------------------------------------------------
    auto x3 = S.LUSolver(false);
    cout << "Solution by LU Decomposition: " << endl;
    cout << x3 << endl << endl;
    //----------------------------------------------------------------
    // Get solution by LU Decomposition Pivoting
    //----------------------------------------------------------------
    auto x4 = S.LUSolver(true);
    cout << "Solution by LU Decomposition Pivoting: " << endl;
    cout << x4 << endl << endl;
    //----------------------------------------------------------------
    // Is A matrix diagonally dominant?
    //print a clear message
    //----------------------------------------------------------------
    cout << "Is A matrix diagonally dominant? " << endl;
    if(FCmatrixAlgo::DiagonalDominant(A))
        cout << "Yes" << endl << endl;
    else
        cout << "No" << endl << endl;

    //----------------------------------------------------------------
    //Get Solution by Gauss Seidel iterative method
    //----------------------------------------------------------------
    vetor_coluna solSeidel(dim);
    int itmax = 100;
    S.IterativeGaussSeidelSolver(solSeidel, itmax);
    cout << "Solution by Gauss Seidel iterative method: " << endl;
    cout << solSeidel << endl << endl;
    //----------------------------------------------------------------
    //Get Solution by Jacobi iterative method
    //----------------------------------------------------------------
    vetor_coluna solJacobi(dim);
    S.IterativeGaussSeidelSolver(solJacobi, itmax);
    cout << "Solution by Jacobi iterative method: " << endl;
    cout << solJacobi << endl << endl;
    //----------------------------------------------------------------
    //Get determinant of coefficient matrix A 
    //----------------------------------------------------------------
    double det = FCmatrixAlgo::Determinante(A);
    cout << "Determinant of coefficient matrix A: " << endl;
    cout << det << endl << endl;
    return 0;
}