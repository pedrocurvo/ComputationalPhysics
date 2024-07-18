#ifndef __FCMATRIXALGO_H__
#define __FCMATRIXALGO_H__
#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
using namespace std;

class FCmatrixAlgo {
    public:
        FCmatrixAlgo() = default; // compiler do it
        ~FCmatrixAlgo() = default;
        /*
        Implements Gauss elimination
        */
       
        static void GaussEliminationMax(
                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>&,
                    Eigen::Matrix<double,Eigen::Dynamic,1>&
                    ); //no pivoting
        static void GaussElimination(
                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>&,
                    Eigen::Matrix<double,Eigen::Dynamic,1>&
                    ); //no pivoting
        static void GaussEliminationPivotAlways(
                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>&,
                    Eigen::Matrix<double,Eigen::Dynamic,1>&,
                    Eigen::Matrix<double,Eigen::Dynamic,1>& // row order indexing
                    ); //make pivoting

        static void GaussEliminationPivot(
                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>&,
                    Eigen::Matrix<double,Eigen::Dynamic,1>&,
                    Eigen::Matrix<double,Eigen::Dynamic,1>& // row order indexing
                    ); //make pivoting
        /*            
        Implements LU decomposition (Doolitle)
        */
        static void LUdecomposition(
                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>&, //matrix coeff
                    Eigen::Matrix<double,Eigen::Dynamic,1>&, // row order indexing
                    bool bpivot=false // activate pivoting
                    );
        /*            
        Implements QR decomposition
        */            
        static void QRdecomposition(
                    Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>&,
                    Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>&,
                    Eigen::Matrix<double, Eigen::Dynamic,Eigen::Dynamic>&
                    );
        
        // Auxiliary Functions -------------------------------------------------
        static double maxrowcoeff(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat, int index);


};
Eigen::Matrix<double,Eigen::Dynamic,1> RowOrderIndexing(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat);

#endif