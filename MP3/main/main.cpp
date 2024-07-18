#include <iostream>
#include <Eigen/Eigen>
#include <EqSolver.h>
#include <FCmatrixAlgo.h>
using namespace std;

int main(){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matr(3, 3);
    
    Eigen::Matrix<double, Eigen::Dynamic, 1> vec(3);
    Eigen::Matrix<double, Eigen::Dynamic, 1> vect(3);
    Eigen::Matrix<int, Eigen::Dynamic, 1> vecto(3);
    
    matr << 1, 2, 3,
            4, 5, 6,
            7, 8, 9;
    vec << 1, 5, 9;
    // Eigen::Matrix<double,Eigen::Dynamic,1> vecking = RowOrderIndexing(matr);
    Eigen::Matrix<double, Eigen::Dynamic, 1> index(3);
    index << 1, 2, 3;

    FCmatrixAlgo sol;
    sol.GaussEliminationMax(matr, vec);
    cout << matr << endl;
    cout << index << endl;
    EqSolver eq(matr, vec);
    cout << eq << endl;
    // //sol.GaussElimination(matr, vec);
    // //sol.GaussEliminationPivot(matr, vec, vecking);
    // //sol.LUdecomposition2(matr, vecking, true); //pivoting enabled is accurate and secure
    // //eq.GaussSolver(true);
    // //eq.LUSolver(false);
    // Eigen::Matrix<double, Eigen::Dynamic, 1> solution(3);
    // solution << 0, 0, 0;
    // //int value = 1000;
    // //eq.IterativeGaussSeidelSolver(solution, value, 1e-100);
    // // eq.IterativeJacobiSolver(solution, value, 1e-100);
    // cout << solution << endl;
    // // cout << matr.inverse() * vec << endl;
    // //cout << matr << endl;
    // auto Q = matr;
    // auto R = matr;
    // FCmatrixAlgo::QRdecomposition(matr, Q, R);
    // cout << Q << endl << endl;
    // cout << R << endl << endl;


    return 0;
}
