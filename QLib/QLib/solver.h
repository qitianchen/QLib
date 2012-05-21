/*!
 * \file solver.h
 *
 * \brief Implementation of non-linear solvers
 *
 * The equation should be of the following format, \n
 * y = f(x)
 *
 *  Created on: May 2, 2012
 *      Author: Qitian Chen
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <cmath>
#include "common.h"

namespace QLib
{
	namespace Solver
	{
		/*!
		 * \brief Bisection solver
		 *
		 * This solver requires that the function is monotonic increasing
		 *
		 * \tparam T The equation
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param y The desired value, left hand side of the equation
		 * \param x1 Lower boundary
		 * \param x2 Upper boundary
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The object contains the equation
		 * \return The root of the equation
		 */
		template<class T, double (T::*f)(double) const>
		double Bisection(double y, double x1, double x2,
				double tol, int maxIterations, const T& obj)
		{
			double x0, y0;
			int count = 0;

			if ( x1 == x2 )
			{
				return x1;
			}
			else if ( x1 > x2 )
			{
				// throw some exception
			}

			if ( maxIterations <= 0 )
			{
				maxIterations = MAX_ITERATIONS;
			}

			do
			{
				x0 = (x1+x2)/2;
				y0 = (obj.*f)(x0);
				if ( y0 > y )
				{
					x2 = x0;
				}
				else
				{
					x1 = x0;
				}

			} while ( std::fabs(y0-y) > tol && ++count < maxIterations );

			return x0;
		}

		/*!
		 * \brief Bisection solver
		 *
		 * This solver requires that the function is monotonic increasing
		 * This solver is the same as the one above but assumes that the object
		 * implements the operator () overloading for the right hand side
		 * of the equation.
		 *
		 * \tparam T The equation object
		 *
		 * \param y The desired value, left hand side of the equation
		 * \param x1 Lower boundary
		 * \param x2 Upper boundary
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param f The object implements the function, right hand side of the equation
		 * \return The root of the equation
		 */
		template<class T>
		double Bisection(double y, double x1, double x2,
				double tol, int maxIterations, const T& f)
		{
			double x0, y0;
			int count = 0;

			if ( x1 == x2 )
			{
				return x1;
			}
			else if ( x1 > x2 )
			{
				// throw some exception
			}

			if ( maxIterations <= 0 )
			{
				maxIterations = MAX_ITERATIONS;
			}

			do
			{
				x0 = (x1+x2)/2;
				y0 = f(x0);
				if ( y0 > y )
				{
					x2 = x0;
				}
				else
				{
					x1 = x0;
				}

			} while ( std::fabs(y0-y) > tol && ++count < maxIterations );

			return x0;
		}


		/*!
		 * \brief Newton-Raphson solver
		 *
		 * This implementation is good for the equations having analytic form of
		 * the derivative.
		 *
		 * \tparam T The equation
		 * \tparam (T::*f) The right hand side of the equation
		 * \tparam (T::*fx) The derivative of the right hand side of the equation to x
		 *
		 * \param y The desired value, left hand side of the equation
		 * \param x0 Initial guess of the root
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The object contains the equation
		 * \return The root of the equation
		 */
		template<class T,
		double (T::*f)(double) const,
		double (T::*fx)(double) const>
		double NewtonRaphson(double y, double x0,
				double tol, int maxIterations, const T& obj)
		{
			double y0 = (obj.*f)(x0);
			int count = 0;

			if ( maxIterations <= 0 )
			{
				maxIterations = MAX_ITERATIONS;
			}

			while ( std::fabs(y0-y) > tol && count++ < maxIterations )
			{
				x0 -= y0/(obj.*fx)(x0);
				y0 = (obj.*f)(x0);
			}

			return x0;
		}


		/*!
		 * \brief Secant solver
		 *
		 * This implementation is suitable for problems which are very difficult
		 * to derivative analytic form of derivative.
		 *
		 * \tparam T The equation
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param y The desired value, left hand side of the equation
		 * \param x0 Initial guess of the root
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The object contains the equation
		 * \return The root of the equation
		 */
		template<class T, double (T::*f)(double) const>
		double Secant(double y, double x0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double y0 = (obj.*f)(x0);
			double fx;
			int count = 0;

			if ( maxIterations <= 0 )
			{
				maxIterations = MAX_ITERATIONS;
			}

			while ( std::fabs(y0-y) > tol && count++ < maxIterations )
			{
				fx = ((obj.*f)(x1) - y0)/(x1 - x0);
				x1 = x0;
				x0 -= y0/fx;
				y0 = (obj.*f)(x0);
			}

			return x0;
		}


		/*!
		 * \brief Secant solver
		 *
		 * This implementation is suitable for problems which are very difficult
		 * to derive analytic form of derivative.
		 *
		 * \tparam T The equation
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param y The desired value, left hand side of the equation
		 * \param x0 Initial guess of the root
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param f The object implements the function, right hand side of the equation
		 * \return The root of the equation
		 */
		template<class T>
		double Secant(double y, double x0, double x1,
				double tol, int maxIterations, const T& f)
		{
			double y0 = f(x0);
			double fx;
			int count = 0;

			if ( maxIterations <= 0 )
			{
				maxIterations = MAX_ITERATIONS;
			}

			while ( std::fabs(y0-y) > tol && count++ < maxIterations )
			{
				fx = (f(x1) - y0)/(x1 - x0);
				x1 = x0;
				x0 -= y0/fx;
				y0 = f(x0);
			}

			return x0;
		}


		/*!
		 * \brief Fixed point iteration solver
		 *
		 * Assuming a<= x <= b and a<=g(x)<=b, then the equation x=g(x) has at
		 * least one solution (fixed point) in the interval [a,b]. Thus the
		 * equation f(x) = 0 could be algebraically manipulated into x=g(x).
		 *
		 * \tparam T The equation
		 * \tparam (T::*f) The right hand side of the equation
		 * \tparam (T::*g) The mapping function
		 *
		 * \param y The desired value, left hand side of the equation
		 * \param x0 Initial guess of the root
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The object contains the equation
		 * \return The root of the equation
		 */
		template<class T,
		double (T::*f)(double) const,
		double (T::*g)(double) const>
		double FixedPointIteration(double y, double x0,
				double tol, int maxIterations, const T& obj)
		{
			double y0 = (obj.*f)(x0);
			int count = 0;

			while ( std::fabs(y0-y) > tol && count++ < maxIterations )
			{
				x0 = y0/(obj.*g)(x0);
				y0 = (obj.*f)(x0);
			}

			return x0;
		}

	}
}

#endif /* SOLVER_H_ */
