#ifndef __FCMATRIXALGO_H__
#define __FCMATRIXALGO_H__

#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
using namespace std;

// typedefs:
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matriz_din;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> vetor_coluna;

class FCmatrixAlgo
{
    public:
        // Constructor: ----------------------------------------------------------------------------

        FCmatrixAlgo() = default; // compiler do it
        ~FCmatrixAlgo() = default;

        // Implements Gauss elimination: -----------------------------------------------------------

        /* 
        Esta função recebe uma matriz quadrada e um vetor e aplica-lhes o método de Gauss dando
        origem a uma matriz triangular superior. Em particular, esta função realiza a troca de
        linhas de modo a que o pivot seja sempre o elemento máximo de cada coluna:
        */
        static void GaussEliminationMax(
                    matriz_din&,
                    vetor_coluna&
                    ); //no pivoting

        /*
        Esta função recebe uma matriz quadrada e um vetor e palica-lhes o método de Gauss originando
        uma matriz triangular superior. Em particular, esta função realiza a troca de linhas apenas
        quando o pivot original de cada linha é nulo, ou seja, não é válido:
        */
        static void GaussElimination(
                    matriz_din&,
                    vetor_coluna&
                    ); //no pivoting
        
        /*
        Esta função recebe uma matriz quadrada e um vetor e aplica-lhes o método de Gauss com 
        pivotagem dando origem a uma matriz triangular superior. Em particular, esta função aplica
        o método da pivotagem para cada uma das iterações sobre as linhas da matriz, verificando 
        qual é o pivot que tem o maior tamanho relativo:
        */
        static void GaussEliminationPivotAlways(
                    matriz_din&,
                    vetor_coluna&,
                    vetor_coluna& // row order indexing
                    ); //make pivoting

        /*
        Esta função recebe uma matriz quadrada e um vetor e aplica-lhes o método de Gauss com 
        pivotagem dando origem a uma matriz triangular superior. Em particular, esta função aplica
        o método de pivotagem com os valores máximos de cada linha da matriz inicial, utilizando-os
        para avaliar qual é o pivot com o maior tamnaho relativo:
        */
        static void GaussEliminationPivot(
                    matriz_din&,
                    vetor_coluna&,
                    vetor_coluna& // row order indexing
                    ); //make pivoting
       
        // Implements LU decomposition (Doolitle): -------------------------------------------------
        
        /*
        Esta função recebe uma matriz e um vetor e aplica-lhes o método de Gauss, preenchendo as
        entradas abaixo da diagonal com os fatores de multiplicação que anularam as mesmas. A
        função tem ainda a opção de aplicar o método de Gauss com ou sem pivotagem:
        */
        static void LUdecomposition(
                    matriz_din&, //matrix coeff
                    vetor_coluna&, // row order indexing
                    bool bpivot = false // activate pivoting
                    );
        
        // Implements QR decomposition: ------------------------------------------------------------

        /*
        Esta função recebe uma matriz que vai decompor em duas matrizes novas, a matriz Q e a matriz
        R, que também são recebidas por referência:
        */       
        static void QRfactorization(
                    matriz_din&,
                    matriz_din&,
                    matriz_din&
                    );

        // Iterative QR: ---------------------------------------------------------------------------

        /*
        Esta função recebe uma matriz, faz a sua decomposição QR, descobre os seus valores próprios
        e preenche uma matriz com estes valores na diagonal:
        */
        static void IterativeQR(
                    matriz_din&,
                    matriz_din&,
                    matriz_din&,
                    double
                    );

        // Determinant: ----------------------------------------------------------------------------

        /*
        Esta função recebe uma matriz e devolve o valor do determinante da matriz:
        */
        static double Determinante(
                    matriz_din&
                    );

        // Inverse matrix: -------------------------------------------------------------------------
        
        /*
        Esta função recebe duas matrizes, uma delas vazia e outra que pretendemos inverter. 
        A função enche a matriz vazia com as entradas da matriz inversa da outra:
        */
        static void Invert(
                    matriz_din,
                    matriz_din&
                    );
        
        /*
        Esta função recebe uma matriz e avalia se a matriz tem ou não dominância diagonal:
        */
        static bool DiagonalDominant(
                    matriz_din&
                    );
        
        // Auxiliary Functions: --------------------------------------------------------------------

        /*
        Esta função recebe uma matriz e o index para uma das linhas da matriz, devolvendo a entrada
        máxima em módulo dessa linha:
        */
        static double maxrowcoeff(
                    matriz_din&,
                    int index
                    );

};

/*
Esta função recebe uma matriz e devolve um vetor que possui a entrada máxima em módulo de cada
uma das linhas da matriz:
*/
vetor_coluna RowOrderIndexing(matriz_din& mat);

#endif