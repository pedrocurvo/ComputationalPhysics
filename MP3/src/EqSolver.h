#ifndef __EQSOLVER_H__
#define __EQSOLVER_H__
#include <Eigen/Dense>
using namespace std;


class EqSolver {
    public:
        // constructors and destructor
        EqSolver() = default;
        EqSolver(
                const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>&, //matrix coeffs
                const Eigen::Matrix<double,Eigen::Dynamic,1>& // vector of constants
            );
        ~EqSolver() = default;
        
        // output
        friend ostream& operator<<(ostream&, const EqSolver&);
        
        // solvers
        const Eigen::Matrix<double,Eigen::Dynamic,1>& GaussSolver(bool pivot=false);
        const Eigen::Matrix<double,Eigen::Dynamic,1>& LUSolver(bool pivot=false);
        const Eigen::Matrix<double,Eigen::Dynamic,1>& MatrixInversion();


        void IterativeJacobiSolver(
            Eigen::Matrix<double,Eigen::Dynamic,1>&, // starting solution
            int& itmax, //nb of max iterations
            double tol=1.E-3); // tolerance on convergence
        void IterativeGaussSeidelSolver(
            Eigen::Matrix<double,Eigen::Dynamic,1>&,
            int& itmax,
            double tol=1.E-3);
            
    private:
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> M; // coefficients matrix
        Eigen::Matrix<double,Eigen::Dynamic,1> b; // constants vector
};

#endif