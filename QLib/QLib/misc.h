/*!
 * \file misc.h
 *
 * \brief Utilities
 *
 *  Created on: May 9, 2012
 *      Author: Administrator
 */

#ifndef MISC_H_
#define MISC_H_

#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace QLib
{
	namespace Misc
	{
		/*!
		 * \brief LU Factorization on a tridiagonal matrix
		 *
		 * Tridiagonal matrix is very special and can be LU decomposed easily.
		 *
		 * \param A The tridiagonal matrix
		 * \param L Lower matrix, with diagonal elements be 1
		 * \param U Upper matrix
		 */
		void TridiagonalLUFactorization(const MatrixXd& A, MatrixXd& L, MatrixXd& U)
		{
			// check whether its tridiagonal
			if ( A.rows() != A.cols() )
			{
				// throw exception
			}

			int N = A.rows();
			int i, j;

			for ( j = 2; j < N; j++ )
			{
				if ( A(0, j) != 0 )
				{
					// throw exception
				}
			}
			for ( i = 1; i < N-1; i++ )
			{
				for ( j = 0; j < i; j++ )
				{
					if ( A(i,j) != 0 )
					{
						// throw exception
					}
				}
				for ( j = i+3; j < N; j++ )
				{
					if ( A(i,j) != 0 )
					{
						// throw exception
					}
				}
			}
			for ( j = 0; j < N-2; j++ )
			{
				if ( A(N-1,j)!= 0 )
				{
					// throw exception
				}
			}

			L = U = MatrixXd::Zero(N,N);
			L(0,0) = 1; U(0,0) = A(0,0); U(0,1) = A(0,1);
			for ( i = 1; i < N-1; i++ )
			{
				L(i,i-1) = A(i,i-1)/U(i-1,i-1);
				L(i,i) = 1;
				U(i,i) = A(i,i)-L(i,i-1)*A(i-1,i);
				U(i,i+1) = A(i,i+1);
			}
			L(N-1, N-2) = A(N-1,N-2)/U(N-2,N-2); L(N-1,N-1) = 1;
			U(N-1, N-1) = A(N-1,N-1)-L(N-1,N-2)*A(N-2,N-1);

			//cout << "LU is " << endl << L*U << endl;
			//cout << "A is " << endl << A << endl;

		}
	}
}

#endif /* MISC_H_ */
