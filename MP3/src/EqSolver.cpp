#include <EqSolver.h>
#include <Eigen/Dense>
#include <FCmatrixAlgo.h>
#include <iostream>
using namespace std;

EqSolver::EqSolver(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat,
                   const Eigen::Matrix<double,Eigen::Dynamic,1>& vec // vector of constants
            )
            {
              M = mat;
              b = vec;
            };


const Eigen::Matrix<double,Eigen::Dynamic,1>& EqSolver::GaussSolver(bool pivot){
    Eigen::Matrix<double,Eigen::Dynamic,1> solution(b.rows());
    FCmatrixAlgo sol;

    if(pivot)
    {
        Eigen::Matrix<double,Eigen::Dynamic,1> vec = RowOrderIndexing(M);
        sol.GaussEliminationPivot(M, b, vec);
    }else{
        sol.GaussElimination(M, b);
    };
    for(int i = M.cols() -1; i >= 0; i--)
    {
        for(int j = i - 1; j >= 0; j--){
            double ratio = M(j, i) / M(i, i);
            b(j) -= ratio * b(i);
            M(j, i) = 0; // no need to subtract, loose precision, just replace by zero
        }
    };
    for(int row = 0; row < M.rows(); ++row)
    { 
        solution(row) = b(row) / M(row, row);
    };
    return solution;
}


const Eigen::Matrix<double,Eigen::Dynamic,1>& EqSolver::LUSolver(bool pivot){
    FCmatrixAlgo sol;
    sol.LUdecomposition(M, b, pivot);
    Eigen::Matrix<double,Eigen::Dynamic,1> x(b.rows());
    Eigen::Matrix<double,Eigen::Dynamic,1> y(b.rows());
    for (int k=0; k<M.rows(); k++)
    {
        double sumC = 0.;
        for (int i=0; i<k; i++) {
            sumC += y(i) * M(k, i);
        }
        y(k) = b(k) - sumC;
    }

    for (int k=M.rows()-1; k>=0; k--)
    {
        double sumC = 0.;
        for (int i=k+1; i< M.rows(); i++) {
            sumC += x(i)*M(k, i);
        }
        x(k) = (y(k) - sumC)/M(k, k);
    }
    return x;
};


const Eigen::Matrix<double,Eigen::Dynamic,1>& EqSolver::MatrixInversion(){
    Eigen::Matrix<double,Eigen::Dynamic,1> x = M.inverse() * b;
    return x;
};


void EqSolver::IterativeJacobiSolver(
    Eigen::Matrix<double, Eigen::Dynamic, 1> &solution , // vector coefficients
    int &N, // number of iterations
    double tol // tolerance
){
    auto R(solution);
    
    auto Raux(R);
    int niter = 0;
    bool good = false;

    while(niter < N && (!good)){
        for(int i = 0; i < M.rows(); i++){
            double sum = 0;
            for(int j = 0; j < M.cols(); j++){
                if(i != j) sum += M(i, j) * R(j, 0);
                
            }
            
            Raux(i, 0) = 1/M(i, i) * (b(i) - sum);
        };

        auto vdif = Raux - R;

        for(int k = 0; k < vdif.rows(); k++){
            if(abs(vdif(k, 0)) < tol){
                good = true;
                }else{
                    good = false;
                }
        };

        R = Raux;
        niter++;
    }

    N = niter;
    solution = R;
}


void EqSolver::IterativeGaussSeidelSolver(
            Eigen::Matrix<double,Eigen::Dynamic,1>& m,
            int& itmax,
            double tol){
    // linear system of m unknowns
   Eigen::Matrix<double,Eigen::Dynamic,1> x(m.rows()); //zero's
   Eigen::Matrix<double,Eigen::Dynamic,1> x_aux(m.rows()); //zero's
    bool btol = false;
    int it = 0.;
    while (!btol && (it++ < 1000)) {
        x_aux = x;
        for (int i=0; i< m.rows(); i++) {
            x[i] = 0.;
            for (int j=0; j< m.rows(); j++){
                if (i != j) x(i) += -M(i, j) * x(j);
            }
            x(i) += b(i);
            x(i) /= M(i, i);
    //guarantee that all vector entries are converging equally
            if (abs(x(i) - x_aux(i)) < tol){
                btol = true;
            }else{
                btol = false;
            }
        }
        
    }
    m = x;
}

ostream& operator<<(ostream&, const EqSolver& E){
    cout << "Matrix: " << endl << E.M << endl << "Vector: " << endl << E.b << endl;
    return cout;
}