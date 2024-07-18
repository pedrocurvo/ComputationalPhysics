#include <FCmatrixAlgo.h>
#include <iostream>
using namespace std;

void FCmatrixAlgo::GaussEliminationMax(matriz_din& mat,vetor_coluna& vector)
{
	for (int k = 0; k < mat.rows() - 1; k++)
	{
		// Initialize maximum value and index for pivot
		int i_max = k;
		double v_max = mat(i_max, k);

		//find the greatest pivot always
		for (int line = k + 1; line < mat.rows(); line++)
		{
			if (abs(mat(line, k)) > v_max) v_max = mat(line, k), i_max = line;
		}

		// If after the code before, the pivot is 0 then matrix is singular
		if (!mat(k,i_max))
			throw logic_error("Matrix is singular");

		//Swap the maximum pivot row with row choosen before
		if (i_max != k)
			mat.row(k).swap(mat.row(i_max));
			vector.row(k).swap(vector.row(i_max));

		for (int line = k + 1; line < mat.rows(); line++)
		{
			//factor of pivot_second_row/pivot_leading_row
			double ratio = mat(line, k) / mat(k, k);

			//subtract the ratio times the kth row from the current row
			for (int col = k + 1; col < mat.cols(); col++)
			{
				mat(line, col) -= mat(k, col) * ratio;
			}
			vector.row(line) -= vector.row(k) * ratio;

			// fills with zeros below the diagonal
			mat(line, k) = 0;
		}
	}
}


void FCmatrixAlgo::GaussElimination(matriz_din& mat, vetor_coluna& vec)
{	
	for(int line = 0; line < mat.rows(); line++)
	{
		if(mat(line, line) == 0){
			for(int col = line + 1; col < mat.cols(); col++)
			{
				if(mat(col, line) != 0)
				{
					mat.row(line).swap(mat.row(col));
					vec.row(line).swap(vec.row(col));
					break;
				}

			}
		}
		if(mat(line, line) == 0) throw logic_error("Matrix is singular");
		for(int col = line + 1; col < mat.cols(); col++)
		{
			double ratio = mat(col, line) / mat(line, line);
			mat.row(col) -= mat.row(line) * ratio;
			vec(col) -= ratio * vec(line);
			mat(col, line) = 0;

		}
	}
}


double FCmatrixAlgo::maxrowcoeff(matriz_din& mat,int index)
{
	double coeff = 0;
	for(int col = 0; col < mat.cols(); col++)
	{
		if(abs(mat(index, col)) > coeff) coeff = abs(mat(index, col));
	}
	return coeff;
}


Eigen::Matrix<double,Eigen::Dynamic,1> RowOrderIndexing(matriz_din& mat)
{
	vetor_coluna vec(mat.rows());
	for(int line = 0; line < mat.rows(); line++)
	{
		vec(line) = FCmatrixAlgo::maxrowcoeff(mat, line);
	}
	return vec;
}


void FCmatrixAlgo::GaussEliminationPivotAlways(matriz_din& mat,vetor_coluna& vector,vetor_coluna& vec_max_coeffs_row)
{
	for (int k = 0; k < mat.rows() - 1; k++)
	{

		// Find the rows to swap:
		std::vector <double> index_controler;
		int z = 0;
		for(int line = k; line < mat.rows(); line++)
		{
			if(mat.row(line).maxCoeff() != 0)
			{
				vec_max_coeffs_row(z) = abs(mat(line, k) / maxrowcoeff(mat, line));
				index_controler.push_back(z + k);
			}
			else
			{
				vec_max_coeffs_row(z) = 0; // to avoid division by zero
				index_controler.push_back(z + k);
			}
			z++;
		}

		int swap_index = k;
		double max_coeff = vec_max_coeffs_row(0);
		for(int i = 0; i < index_controler.size(); i++)
		{
			if(vec_max_coeffs_row(i) > max_coeff)
			{
				max_coeff = vec_max_coeffs_row(i);
				swap_index = index_controler[i];
			}
		}

		// Check if its singular
		int i_max = k;
		double v_max = mat(i_max, k);

		//find the greatest pivot always
		for (int line = k + 1; line < mat.rows(); line++)
		{
			if (abs(mat(line, k)) > v_max) v_max = mat(line, k), i_max = line;
		}

		// If after the code before, the pivot is 0 then matrix is singular
		if (!mat(k,i_max)) throw logic_error("Matrix is singular");

		//Swap the maximum pivot row with row choosen before
		if (swap_index != k)
		{
			mat.row(k).swap(mat.row(swap_index));
			vector.row(k).swap(vector.row(swap_index));
		}

		for (int line = k+1; line < mat.rows(); line++)
		{
			//factor of pivot_second_row/pivot_leading_row
			double ratio = mat(line, k) / mat(k, k);

			//subtract the ratio times the kth row from the current row
			for (int col = k+1; col < mat.rows(); col++)
			{
				mat(line, col) -= mat(k , col) * ratio;
			}
			vector.row(line) -= vector.row(k) * ratio;

			// fills with zeros below the diagonal
			mat(line, k) = 0;
		}

	vec_max_coeffs_row.resize(vec_max_coeffs_row.rows() - 1, 1);
	index_controler.clear();
	}
}


void FCmatrixAlgo::GaussEliminationPivot(matriz_din& mat,vetor_coluna& vector,vetor_coluna& row_order)
{
	vetor_coluna vec_mxcoeffs = RowOrderIndexing(mat);
	for(int i = 0; i < vec_mxcoeffs.rows(); i++)
	{
		if(vec_mxcoeffs(i) == 0)
		{
			throw logic_error("Matrix is singular");
		}
	}
	for(int row = 0; row < mat.rows(); ++row)
	{
		double pivot_control = mat(row, row) / vec_mxcoeffs(row);
		int index_controler = row;
		for(int k = row; k < mat.rows(); k++)
		{
			if( (abs(mat(k, row) / vec_mxcoeffs(k))) > pivot_control)
			{
				pivot_control = abs(mat(k, row) / vec_mxcoeffs(k));
				index_controler = k;
			}
		}
		// Swaps:
		if(index_controler != row)
		{
			mat.row(row).swap(mat.row(index_controler));
			vector.row(row).swap(vector.row(index_controler));
			row_order.row(row).swap(row_order.row(index_controler));
		}

		// Check if its singular
		int i_max = row;
		double v_max = mat(i_max, row);
		//find the greatest pivot always
		for (int line = row + 1; line < mat.rows(); line++)
		{
			if (abs(mat(line, row)) > v_max) v_max = mat(line, row), i_max = line;
		}
		// If after the code before, the pivot is 0 then matrix is singular
		if (!mat(row,i_max)) throw logic_error("Matrix is singular");
		
		for (int line = row + 1; line < mat.rows(); line++)
		{
			//factor of pivot_second_row/pivot_leading_row
			double ratio = mat(line, row) / mat(row, row);

			//subtract the ratio times the kth row from the current row:
			mat.row(line) -= mat.row(row) * ratio;
			vector.row(line) -= vector.row(row) * ratio;

			// fills with zeros below the diagonal, avoid loosing precision by subtracting
			mat(line, row) = 0;
		}
	}
}


void FCmatrixAlgo::LUdecomposition(matriz_din& mat,vetor_coluna& vec,bool bpivot)
{
	vetor_coluna vec_mxcoeffs = RowOrderIndexing(mat);
	if(bpivot)
	{
		
		for(int i = 0; i < vec_mxcoeffs.rows(); i++)
		{
			if(vec_mxcoeffs(i) == 0){
				throw logic_error("Matrix is singular");
			}
		}
	}

	for(int row = 0; row < mat.rows(); ++row)
	{
		if(bpivot)
		{
			double pivot_control = mat(row, row) / vec_mxcoeffs(row) ;
			int index_controler = row;
			for(int k = row; k < mat.rows(); k++)
			{
				if( (abs(mat(k, row) / vec_mxcoeffs(k))) > pivot_control)
				{
					pivot_control = abs(mat(k, row) / vec_mxcoeffs(k));
					index_controler = k;
				}
			}
			
			if(index_controler != row){
				mat.row(row).swap(mat.row(index_controler));
				vec.row(row).swap(vec.row(index_controler));
			}
		}

		// Check if its singular
		int i_max = row;
		double v_max = mat(i_max, row);
		//find the greatest pivot always
		for (int i = row + 1; i < mat.rows(); i++)
		{
			if (abs(mat(i, row)) > v_max) v_max = mat(i, row), i_max = i;
		}
		// If after the code before, the pivot is 0 then matrix is singular
		if (!mat(row,i_max)) throw logic_error("Matrix is singular");

		
		for (int i = row + 1; i < mat.rows(); i++)
		{
			//factor of pivot_second_row/pivot_leading_row
			double ratio = mat(i, row) / mat(row, row);
			//subtract the ratio times the kth row from the current row
			mat.row(i) -= mat.row(row) * ratio;
			// fills the lower diagonal part with LOWER matrix coefficients
			mat(i, row) = ratio;
		}
	}
}


void FCmatrixAlgo::QRfactorization(matriz_din& mat,matriz_din& Q,matriz_din& R)
{
		Eigen::MatrixXd Q_aux = mat;
		for(int line = 1; line < mat.rows(); line++)
		{
			for(int col = line - 1; col >= 0; col--)
			{
				Q_aux.col(line) -= ((Q_aux.col(col).transpose() * Q_aux.col(line))/(Q_aux.col(col).transpose() * Q_aux.col(col)))(0) * Q_aux.col(col);
			}
			Q_aux.col(line).normalize();

		}
		Q_aux.col(0).normalize();
		R = Q_aux.transpose() * mat;
		for(int col = 0; col < R.cols(); col++)
		{
			for(int line = col + 1; line < R.rows(); line++)
			{
				R(line,col) = 0;
			}
		}
		Q = Q_aux;
}


bool FCmatrixAlgo::DiagonalDominant(matriz_din& mat)
{
	bool res = false;
	for(int line = 0; line < mat.rows(); line++)
	{	
		double sum = 0;
		for(int col = 0; col < mat.cols(); col++){
			sum += abs(mat(line,col));
		}
		if(sum > mat(line, line)) res = true;
	}
	return res;
}


double FCmatrixAlgo::Determinante(matriz_din& mat)
{	
	Eigen::Matrix<double,Eigen::Dynamic, 1> vector(mat.rows());
    GaussElimination(mat, vector);
    double determinante = 1;
    for(int line = 0; line < mat.rows(); ++line)
	{
        determinante *= mat(line, line);
    }
	return determinante;
}


void FCmatrixAlgo::Invert(matriz_din mat,matriz_din& mat_inverted)
{
	// Vamos transformar a matriz inversa na matriz indentidade com o tmanaho da matriz que 
	// queremos inverter:

	mat_inverted.resize(mat.rows(), mat.cols());
	for (int row = 0; row < mat_inverted.rows(); ++row)
	{
		mat_inverted(row, row) = 1;
	}
	
	// Vamos aplicar o metodo de Gauss Jordan a ambas as matrizes simultaneamente:

	// Gauss:

    for (int i = 0; i < mat.rows(); ++i)
	{
        double pivot = mat(i,i);

        // só trocamos linhas quando o pivot é igual a 0:
        if (mat(i,i) == 0)
		{
            std::vector <double> coluna;
            for (int col = i; col < mat.cols(); ++col) coluna.push_back(abs(mat(col,i)));
            pivot = *max_element(coluna.begin(), coluna.end()); //pivot is the maximum of each column
            if (pivot == 0) throw logic_error ("singular matrix: no inverse");

			// swap lines:
			for (int line = i + 1; line < mat.rows(); ++line)
			{
				if (abs(mat(line,i)) == pivot && mat(line,i) != 0)
				{
					mat.row(i).swap(mat.row(line)); //swap when pivot is found
					mat_inverted.row(i).swap(mat_inverted.row(line));
				}           
			}
        }

        /* Vamos percorrer as linhas abaixo da linha iterada e vamos subtrair-lhe um multiplo da 
		linha com pivot de modo a que essa coluna seja nula a partir do pivot: */
        for (int line = i + 1; line < mat.rows(); ++line)
		{
            double fator_mul = mat(line,i) / mat(i,i); //find factor for subtraction
			mat.row(line) -= mat.row(i) * fator_mul; //subtract
			mat_inverted.row(line) -= mat_inverted.row(i) * fator_mul;
        }
	}

	// Jordan:

	// Vamos percorrer todas as colunas, começando pela última:
	for (int col = mat.cols() - 1; col >= 0; col--){
		// Vamos percorrer as linhas acima do pivot e subtrair o múltiplo adequado do pivot:
		for (int line = col - 1; line >= 0; line--){
			double fator_mul = mat(line,col) / mat(col,col);
			// mat(line,col) -= mat(col,col) * fator_mul;
			mat.row(line) -= mat.row(col) * fator_mul;
			mat_inverted.row(line) -= mat_inverted.row(col) * fator_mul;
		}
	}

	// Vamos transformar a matriz original na identidade e obter a inversa:
	for (int row = 0; row < mat.rows(); row++){
		double fator_div = mat(row,row);
		mat(row,row) /= fator_div;
		mat_inverted.row(row) /= fator_div;
	}
}


void FCmatrixAlgo::IterativeQR(matriz_din& V,matriz_din& E,matriz_din& X,double eps)
{
	//E fica com os eigenvalues 
	//X fica com os eigenvectors
	matriz_din Q;
	matriz_din R;

	//obter a primeira decomposicao de V
	FCmatrixAlgo::QRfactorization(V,Q,R);
	X = Q;
	bool continua = true;

	while (continua)
	{
		//definir a nova matriz (V1 = R1*Q1) usando o R e o Q que veem da iteracao anterior 
		E = R*Q;

		//fatorizar a nova matriz para criar os novos R e Q
		FCmatrixAlgo::QRfactorization(E,Q,R);

		X *= Q;

		continua = false;
		//check if eps is met
		for (int line = 0; line < E.rows(); line++)
		{
			for (int col = 0; col < E.cols(); col++)
			{
				if (line != col)
				{
					if (abs(E(line,col)) > eps)
					{
						continua = true;
					}
				}
			}
		}
	}
}