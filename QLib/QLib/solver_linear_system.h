/*!
 * \file solver_linear_system.h
 *
 * \brief Linear system solver
 *
 * The file specifies the linear system as follows,\n
 * Ax = b
 *
 *  Created on: May 6, 2012
 *      Author: Qitian Chen
 */

#ifndef SOLVER_LINEAR_SYSTEM_H_
#define SOLVER_LINEAR_SYSTEM_H_

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

namespace QLib
{
	namespace Solver
	{
		/*!
		 * \brief Solve by inverting the left hand size matrix
		 *
		 * \param A The matrix on the left hand side
		 * \param b The vector on the right hand side
		 * \return x, the solution to the system
		 */
		VectorXd DirectInversion(const MatrixXd& A, const VectorXd& b)
		{
			if ( A.rows() != A.cols() || A.rows() != b.rows() )
			{
				// throw exception
			}

			// check rank?

			 return A.inverse() * b;
		}

		/*!
		 * \brief Solve by Gaussian elimination
		 *
		 * \param A The matrix on the left hand side
		 * \param b The vector on the right hand side
		 * \return x, the solution to the system
		 */
		VectorXd GaussianElimination(const MatrixXd& A, const VectorXd& b)
		{
			if ( A.rows() != A.cols() || A.rows() != b.rows() )
			{
				// throw exception
			}

			// check rank?

			int ncols = A.cols();
			int nrows = A.rows();
			int i, j;

			MatrixXd Ab(nrows, ncols+1);
			VectorXd x = b;

			// nasty copy
			for ( i = 0; i < nrows; i++ )
			{
				for ( j = 0; j < ncols; j++ )
				{
					Ab(i,j) = A(i,j);
				}
				Ab(i, ncols) = b[i];
			}

			int icol;
			double mul, val;
			for ( icol = 0; icol < ncols; icol++ )
			{
				//cout << "Ref row " << icol << " is " << Ab.row(icol) << endl;
				mul = Ab(icol, icol);
				for ( i = icol; i <= ncols; i++ )
				{
					Ab(icol, i) = Ab(icol, i)/mul;
				}

				//cout << "Unitized ref row " << icol << " is " << endl << Ab.row(icol) << endl;

				for ( i = icol+1; i < nrows; i++ )
				{
					//cout << "Row " << i << " is " << endl << Ab.row(i) << endl;
					mul = Ab(i,icol);
					//cout << "Multiplier is " << mul << endl;
					for ( j = icol; j <= ncols; j++ )
					{
						Ab(i,j) -= Ab(icol, j)*mul;
					}
					//cout << "Reduced row " << i << " is " << endl << Ab.row(i) << endl;
					//mul *= 9;
				}
				//cout << "Step "<<icol<< "Ab is "<<endl<<Ab<<endl;
			}

			for ( icol = ncols-1; icol >= 0; icol-- )
			{
				val = 0;
				for ( i = icol+1; i < ncols; i++ )
				{
					val += Ab(icol,i)*x[i];
				}
				x[icol] = Ab(icol, ncols) - val;
			}

			return x;
		}
		

		/*!
		 * \brief Solve by Jacobi iterative method
		 *
		 * \param A The matrix on the left hand side
		 * \param b The vector on the right hand side
		 * \param x0 Initial guess of the solution
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \return x, the solution to the system
		 */
		VectorXd JacobiIterative(const MatrixXd& A, const VectorXd& b,
				const VectorXd& x0, double tol, int maxIterations)
		{
			if ( A.rows() != A.cols() || A.rows() != b.rows() )
			{
				// throw exception
			}

			int count = 0;

			if ( maxIterations <= 0 )
			{
				maxIterations = MAX_ITERATIONS;
			}

			// check rank?
			int nrows = b.rows();
			int ncols = A.cols();
			double dist;

			VectorXd x = x0, x_ = x0, tmpV = x0;
			MatrixXd N = MatrixXd::Zero(ncols, nrows);
			int i;

			for ( i = 0; i < ncols; i++ )
			{
				N(i,i) = A(i,i);
			}

			MatrixXd P = N - A;
			MatrixXd N_ = N.inverse();

			//cout << "A is " << endl << A << endl;
			//cout << "N is " << endl << N << endl;
			//cout << "Inverse of N is " << endl << N_ << endl;
			//cout << "P is " << endl << P << endl;

			do
			{
				//cout << "Before: x is " << endl << x << endl;

				x = N_*(b+P*x);
				tmpV = x - x_;
				dist = 0;
				for ( i = 0; i < ncols; i++ )
				{
					dist += tmpV[i]*tmpV[i];
				}
				dist = sqrt(dist);
				x_ = x;

				//cout << "After: x is " << endl << x << endl;
				//cout << "Diff is " << endl << tmpV << endl;
				//cout << "Error is " << dist << endl;

			} while ( dist > tol && ++count < maxIterations );

			return x;
		}

		/*!
		 * \brief Solve by Gauss-Seibal iterative method
		 *
		 * \param A The matrix on the left hand side
		 * \param b The vector on the right hand side
		 * \param x0 Initial guess of the solution
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \return x, the solution to the system
		 */
		VectorXd GaussSeidal(const MatrixXd& A, const VectorXd& b,
				const VectorXd& x0, double tol, int maxIterations)
		{
			if ( A.rows() != A.cols() || A.rows() != b.rows() )
			{
				// throw exception
			}

			int count = 0;

			if ( maxIterations <= 0 )
			{
				maxIterations = MAX_ITERATIONS;
			}

			// check rank?
			int nrows = b.rows();
			int ncols = A.cols();
			double dist;

			VectorXd x = x0, x_ = x0, tmpV = x0;
			MatrixXd N = MatrixXd::Zero(ncols, nrows);
			int i, j;

			for ( i = 0; i < ncols; i++ )
			{
				for ( j = 0; j <= i; j++ )
				{
					N(i,j) = A(i,j);
				}
			}

			MatrixXd P = N - A;
			MatrixXd N_ = N.inverse();

			//cout << "A is " << endl << A << endl;
			//cout << "N is " << endl << N << endl;
			//cout << "Inverse of N is " << endl << N_ << endl;
			//cout << "P is " << endl << P << endl;

			do
			{
				//cout << "Before: x is " << endl << x << endl;

				x = N_*(b+P*x);
				tmpV = x - x_;
				dist = 0;
				for ( i = 0; i < ncols; i++ )
				{
					dist += tmpV[i]*tmpV[i];
				}
				dist = sqrt(dist);
				x_ = x;

				//cout << "After: x is " << endl << x << endl;
				//cout << "Diff is " << endl << tmpV << endl;
				//cout << "Error is " << dist << endl;

			} while ( dist > tol && ++count < maxIterations );

			return x;
		}

		/*!
		 * \brief Solve by LU factorization
		 *
		 * \param A The matrix on the left hand side
		 * \param b The vector on the right hand side
		 * \return x, the solution to the system
		 */
		VectorXd LUFactorization(const MatrixXd& A, const VectorXd& b)
		{
			if ( A.rows() != A.cols() || A.rows() != b.rows() )
			{
				// throw exception
			}

			return A.lu().solve(b);
		}
	}
}

#endif /* SOLVER_LINEAR_SYSTEM_H_ */
