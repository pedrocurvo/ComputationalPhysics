#include <EqSolver.h>
#include <Eigen/Dense>
#include <FCmatrixAlgo.h>
#include <iostream>
using namespace std;

EqSolver::EqSolver(const matriz_din& mat, const vetor_coluna& vec)
{
    M = mat;
    b = vec;
}


const vetor_coluna EqSolver::GaussSolver(bool pivot)
{
    vetor_coluna solution(b.rows());
    FCmatrixAlgo sol;

    if(pivot)
    {
        vetor_coluna vec = RowOrderIndexing(M);
        sol.GaussEliminationPivot(M, b, vec);
    }
    else
    {
        sol.GaussElimination(M, b);
    }
    for(int col = M.cols() -1; col >= 0; col--)
    {
        for(int line = col - 1; line >= 0; line--){
            double ratio = M(line, col) / M(col, col);
            b(line) -= ratio * b(col);
            M(line, col) = 0; // no need to subtract, loose precision, just replace by zero
        }
    }
    for(int row = 0; row < M.rows(); ++row)
    { 
        solution(row) = b(row) / M(row, row);
    }
    return solution;
}


const vetor_coluna EqSolver::LUSolver(bool pivot)
{
    FCmatrixAlgo sol;
    sol.LUdecomposition(M, b, pivot);
    vetor_coluna x(b.rows());
    vetor_coluna y(b.rows());
    for (int k=0; k<M.rows(); k++)
    {
        double sumC = 0.;
        for (int i = 0; i < k; i++)
        {
            sumC += y(i) * M(k, i);
        }
        y(k) = b(k) - sumC;
    }

    for (int k = M.rows() - 1; k >= 0; k--)
    {
        double sumC = 0.;
        for (int i = k + 1; i < M.rows(); i++)
        {
            sumC += x(i) * M(k, i);
        }
        x(k) = (y(k) - sumC) / M(k, k);
    }
    return x;
}


void EqSolver::IterativeJacobiSolver(vetor_coluna &solution , int &N, double tol)
{
    auto R(solution);
    
    auto Raux(R);
    int niter = 0;
    bool good = false;

    while(niter < N && (!good))
    {
        for(int i = 0; i < M.rows(); i++)
        {
            double sum = 0;
            for(int j = 0; j < M.cols(); j++)
            {
                if(i != j) sum += M(i, j) * R(j, 0);
                
            }
            Raux(i, 0) = 1/M(i, i) * (b(i) - sum);
        }

        auto vdif = Raux - R;

        for(int k = 0; k < vdif.rows(); k++)
        {
            if(abs(vdif(k, 0)) < tol)
            {
                good = true;
                }
                else
                {
                    good = false;
                }
        }

        R = Raux;
        niter++;
    }

    N = niter;
    solution = R;
}


void EqSolver::IterativeGaussSeidelSolver(vetor_coluna & m, int& itmax, double tol)
{
    // linear system of m unknowns
    vetor_coluna x(m.rows()); //zero's
    vetor_coluna x_aux(m.rows()); //zero's
    bool btol = false;
    int it = 0.;
    while (!btol && (it++ < 1000))
    {
        x_aux = x;
        for (int line = 0; line < m.rows(); line++)
        {
            x[line] = 0.;
            for (int col = 0; col < m.rows(); col++)
            {
                if (line != col) x(line) += -M(line, col) * x(col);
            }
            x(line) += b(line);
            x(line) /= M(line, line);
            
            //guarantee that all vector entries are converging equally
            if (abs(x(line) - x_aux(line)) < tol)
            {
                btol = true;
            }
            else
            {
                btol = false;
            }
        }
        
    }
    m = x;
}


ostream& operator<<(ostream&, const EqSolver& E)
{
    cout << "Matrix: " << endl << E.M << endl << "Vector: " << endl << E.b << endl;
    return cout;
}