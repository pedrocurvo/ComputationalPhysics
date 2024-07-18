#include <FCmatrixAlgo.h>
#include <iostream>
using namespace std;



void FCmatrixAlgo::GaussEliminationMax(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat, Eigen::Matrix<double, Eigen::Dynamic, 1>& vector){
	
	for (int k = 0; k < mat.rows() - 1; k++){
	// Initialize maximum value and index for pivot
		int i_max = k;
		double v_max = mat(i_max, k);

		//find the greatest pivot always
		for (int i = k + 1; i < mat.rows(); i++){
			if (abs(mat(i, k)) > v_max) v_max = mat(i, k), i_max = i;
		}

		// If after the code before, the pivot is 0 then matrix is singular
		if (!mat(k,i_max))
			throw logic_error("Matrix is singular");

		//Swap the maximum pivot row with row choosen before
		if (i_max != k)
			mat.row(k).swap(mat.row(i_max));
			vector.row(k).swap(vector.row(i_max));

		for (int i = k+1; i < mat.rows(); i++){
			//factor of pivot_second_row/pivot_leading_row
			double ratio = mat(i, k) / mat(k, k);

			//subtract the ratio times the kth row from the current row
			for (int j = k+1; j < mat.rows(); j++){
				mat(i, j) -= mat(k , j) * ratio;
			}
			vector.row(i) -= vector.row(k) * ratio;

			// fills with zeros below the diagonal
			mat(i, k) = 0;
		}
	}
};


void FCmatrixAlgo::GaussElimination(
                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat,
                    Eigen::Matrix<double,Eigen::Dynamic,1>& vec
                    )
{	
	for(int i = 0; i < mat.rows(); i++)
	{
		if(mat(i, i) == 0){
			for(int j = i + 1; j < mat.cols(); j++)
			{
				if(mat(j, i) != 0)
				{
					mat.row(i).swap(mat.row(j));
					vec.row(i).swap(vec.row(j));
					break;
				}

			}
		}
		if(mat(i, i) == 0) throw logic_error("Matrix is singular");
		for(int j = i + 1; j < mat.cols(); j++)
		{
			double ratio = mat(j, i) / mat(i, i);
			mat.row(j) -= mat.row(i) * ratio;
			vec(j) -= ratio * vec(i);
			mat(j, i) = 0;

		}
	}
};


double FCmatrixAlgo::maxrowcoeff(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat, int index){
		double coeff = 0;
		for(int j = 0; j < mat.cols(); j++){
			if(abs(mat(index, j)) > coeff) coeff = abs(mat(index, j));
		}
		return coeff;
};


Eigen::Matrix<double,Eigen::Dynamic,1> RowOrderIndexing(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat){
	Eigen::Matrix<double,Eigen::Dynamic,1> vec(mat.rows());
	for(int i = 0; i < mat.rows(); i++){
		vec(i) = FCmatrixAlgo::maxrowcoeff(mat, i);
	}
	return vec;
};


void FCmatrixAlgo::GaussEliminationPivotAlways(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat,
                    Eigen::Matrix<double,Eigen::Dynamic,1>& vector,
                    Eigen::Matrix<double,Eigen::Dynamic,1>& vec_max_coeffs_row){
				

	for (int k = 0; k < mat.rows() - 1; k++){

	// Find the rows to swap
		std::vector<double> index_controler;
		int z = 0;
		for(int i = k; i < mat.rows(); i++){
			if(mat.row(i).maxCoeff() != 0){
				vec_max_coeffs_row(z) = abs(mat(i, k) / maxrowcoeff(mat, i));
				index_controler.push_back(z + k);
			}else{
				vec_max_coeffs_row(z) = 0; // to avoid division by zero
				index_controler.push_back(z + k);
			}
			z++;
		};

		int swap_index = k;
		double max_coeff = vec_max_coeffs_row(0);
		for(int i = 0; i < index_controler.size(); i++){
			if(vec_max_coeffs_row(i) > max_coeff){
				max_coeff = vec_max_coeffs_row(i);
				swap_index = index_controler[i];
			}
		}




	// Check if its singular
		int i_max = k;
		double v_max = mat(i_max, k);
		//find the greatest pivot always
		for (int i = k + 1; i < mat.rows(); i++){
			if (abs(mat(i, k)) > v_max) v_max = mat(i, k), i_max = i;
		}
		// If after the code before, the pivot is 0 then matrix is singular
		if (!mat(k,i_max))
			throw logic_error("Matrix is singular");


	//Swap the maximum pivot row with row choosen before
		if (swap_index != k)
			mat.row(k).swap(mat.row(swap_index));
			vector.row(k).swap(vector.row(swap_index));

		for (int i = k+1; i < mat.rows(); i++){
			//factor of pivot_second_row/pivot_leading_row
			double ratio = mat(i, k) / mat(k, k);

			//subtract the ratio times the kth row from the current row
			for (int j = k+1; j < mat.rows(); j++){
				mat(i, j) -= mat(k , j) * ratio;
			}
			vector.row(i) -= vector.row(k) * ratio;

			// fills with zeros below the diagonal
			mat(i, k) = 0;
		}

	vec_max_coeffs_row.resize(vec_max_coeffs_row.rows() - 1, 1);
	index_controler.clear();
	}
};


void FCmatrixAlgo::GaussEliminationPivot(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat,
                    Eigen::Matrix<double,Eigen::Dynamic,1>& vector,
                    Eigen::Matrix<double,Eigen::Dynamic,1>& row_order){
	Eigen::Matrix<double,Eigen::Dynamic,1> vec_mxcoeffs = RowOrderIndexing(mat);
	for(int i = 0; i < vec_mxcoeffs.rows(); i++)
	{
		if(vec_mxcoeffs(i) == 0){
			throw logic_error("Matrix is singular");
		}
	}
	for(int row = 0; row < mat.rows(); ++row)
	{
		
		double pivot_control = mat(row, row) / vec_mxcoeffs(row) ;
		int index_controler = row;
		for(int k = row; k < mat.rows(); k++){
			if( (abs(mat(k, row) / vec_mxcoeffs(k))) > pivot_control){
				pivot_control = abs(mat(k, row) / vec_mxcoeffs(k));
				index_controler = k;
			}
		}
		// Swaps
		if(index_controler != row){
			mat.row(row).swap(mat.row(index_controler));
			vector.row(row).swap(vector.row(index_controler));
			row_order.row(row).swap(row_order.row(index_controler));
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
			// for (int j = row + 1; j < mat.rows(); j++){
			// 	mat(i, j) -= mat(row , j) * ratio;
			// }
			mat.row(i) -= mat.row(row) * ratio;
			vector.row(i) -= vector.row(row) * ratio;

			// fills with zeros below the diagonal, avoid loosing precision by subtracting
			mat(i, row) = 0;
		}
	}
};


void FCmatrixAlgo::LUdecomposition( Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat, //matrix coeff
                    Eigen::Matrix<double,Eigen::Dynamic,1>& vec, // row order indexing
                    bool bpivot) // activate pivoting)
					{
	Eigen::Matrix<double,Eigen::Dynamic,1> vec_mxcoeffs = RowOrderIndexing(mat);
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
		if(bpivot){
			double pivot_control = mat(row, row) / vec_mxcoeffs(row) ;
			int index_controler = row;
			for(int k = row; k < mat.rows(); k++){
				if( (abs(mat(k, row) / vec_mxcoeffs(k))) > pivot_control){
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
};


void FCmatrixAlgo::QRdecomposition(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat,
									Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& Q, 
									Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& R){
		Eigen::MatrixXd Q_aux = mat;
		for(int j = 1; j < mat.rows(); j++){
			for(int i = j - 1; i >= 0; i--){
				Q_aux.col(j) -= ((Q_aux.col(i).transpose() * Q_aux.col(j))/(Q_aux.col(i).transpose() * Q_aux.col(i)))(0) * Q_aux.col(i);
			}
			Q_aux.col(j).normalize();

		}
		Q_aux.col(0).normalize();
		R = Q_aux.transpose() * mat;
		for(int i = 0; i < R.cols(); i++){
			for(int j = i + 1; j < R.rows(); j++){
				R(j,i) = 0;
			}
		}
		Q = Q_aux;
}
