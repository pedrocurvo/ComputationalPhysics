// void FCmatrixAlgo::LUdecomposition(
//                     Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& mat,
//                     Eigen::Matrix<double,Eigen::Dynamic,1>& vector,
//                     bool bpivot
//                     ){

// 	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> lower(mat.rows(), mat.cols());
// 	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> upper(mat.rows(), mat.cols());

// 	for (int i = 0; i < mat.rows(); i++)	{
//         // Upper Triangular
//         for (int k = i; k < mat.rows(); k++){
//             // Summation of L(i, j) * U(j, k)
//             int sum = 0;
//             for (int j = 0; j < i; j++){
//                 sum += (lower(i, j) * upper(j , k));
// 			}
//             // Evaluating U(i, k)
//             upper(i, k) = mat(i, k) - sum;
//         }

// 		// Lower Triangular
// 		for (int k = i; k < mat.rows(); k++)
//         {
//             if (i == k) lower(i, i) = 1; // Diagonal as 1
//             else
//             {
//                 // Summation of L(k, j) * U(j, i)
//                 int sum = 0;
//                 for (int j = 0; j < i; j++)
//                     sum += (lower(k, j) * upper(j, i));
 
//                 // Evaluating L(k, i)
//                 lower(k, i) = (mat(k, i) - sum) / upper(i, i);
//             }
//         }
// 	}

	
        
// 	cout << endl << lower  << endl << endl << upper << endl << endl;
// 	cout << lower * upper << endl;
// };

