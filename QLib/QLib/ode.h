/*!
 * \file ode.h
 *
 * \brief Ordinary differential equation (ODE) solvers
 *
 * The ODE equation should be in the following format \n
 * dy/dx = f(x,y)
 *
 *  Created on: May 2, 2012
 *      Author: Qitian Chen
 */

#ifndef ODE_H_
#define ODE_H_

#include <cmath>
#include "common.h"

namespace QLib
{
	namespace ODE
	{

		/*!
		 * \brief Solve an ODE equation with forward/backward difference method
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T, double(T::*f)(double, double) const>
		double Difference(double x0, double y0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double x0_, y0_, y1, y1_;
			double h;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			if (x0 == x1)
			{
				return y0;
			}

			y1 = y0;
			h = (x1 - x0) / 10;

			x0_ = x0;
			y0_ = y0;

			do
			{
				y1_ = y1;
				while ((x1 - x0) * h > PRECISION)
				{
					y0 += h * (obj.*f)(x0, y0);
					x0 += h;
				}
				y1 = y0;
				h /= 2;
				x0 = x0_;
				y0 = y0_;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

		/*!
		 * \brief Solve an ODE equation with Richardson extrapolation method
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T, double(T::*f)(double, double) const>
		double RichardsonExtrapolation(double x0, double y0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double y1, y1_, x0_h, x0_2h, y0_h, y0_2h;
			double h, h_2;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			if (x0 == x1)
			{
				return y0;
			}

			y1 = y0;

			h = (x1 - x0) / 10;
			h_2 = 2 * h;

			do
			{
				x0_h = x0;
				x0_2h = x0;
				y0_h = y0;
				y0_2h = y0;
				y1_ = y1;
				while ((x1 - x0_2h) * h > PRECISION)
				{
					y0_h += h * (obj.*f)(x0_h, y0_h);
					x0_h += h;
					y0_h += h * (obj.*f)(x0_h, y0_h);
					x0_h += h;
					y0_2h += h_2 * (obj.*f)(x0_2h, y0_2h);
					x0_2h += h_2;
				}
				y1 = 2 * y0_h - y0_2h;
				h /= 2;
				h_2 /= 2;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

		/*!
		 * \brief Solve an ODE equation with implicit method
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T, double(T::*solve)(double, double, double) const>
		double Implicit(double x0, double y0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double x0_, y0_, y1, y1_;
			double h;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			if (x0 == x1)
			{
				return y0;
			}

			y1 = y0;
			h = (x1 - x0) / 10;

			x0_ = x0;
			y0_ = y0;

			do
			{
				y1_ = y1;
				while ((x1 - x0) * h > PRECISION)
				{
					y0 = (obj.*solve)(x0, y0, h);
					x0 += h;
				}
				y1 = y0;
				h /= 2;
				x0 = x0_;
				y0 = y0_;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

		/*!
		 * \brief Solve an ODE equation with Predictor Corrector method
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T, double(T::*f)(double, double) const>
		double PredictorCorrector(double x0, double y0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double x0_, y0_, y0p, y1, y1_;
			double h;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			if (x0 == x1)
			{
				return y0;
			}

			y1 = y0;
			h = (x1 - x0) / 10;

			x0_ = x0;
			y0_ = y0;

			do
			{
				y1_ = y1;
				while ((x1 - x0) * h > PRECISION)
				{
					y0p = y0 + h * (obj.*f)(x0, y0);
					x0 += h;
					y0 += h * (obj.*f)(x0, y0p);
				}
				y1 = y0;
				h /= 2;
				x0 = x0_;
				y0 = y0_;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

		/*!
		 * \brief Solve an ODE equation with trapezium method
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T, double(T::*f)(double, double) const>
		double Trapezium(double x0, double y0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double x0_, y0_, x0p, y0p, y1, y1_, delta_y;
			double h;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			if (x0 == x1)
			{
				return y0;
			}

			y1 = y0;
			h = (x1 - x0) / 10;

			x0_ = x0;
			y0_ = y0;

			do
			{
				y1_ = y1;
				while ((x1 - x0) * h > PRECISION)
				{
					x0p = x0 + h;
					delta_y = (obj.*f)(x0, y0);
					y0p = y0 + h * delta_y;
					y0 += 0.5 * h * (delta_y + (obj.*f)(x0p, y0p));
					x0 = x0p;
				}
				y1 = y0;
				h /= 2;
				x0 = x0_;
				y0 = y0_;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

		/*!
		 * \brief Solve an ODE equation with higher order Taylor series expansion
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 * \tparam (T::*fx) The derivative of right hand side of the equation against x
		 * \tparam (T::*fy) The derivative of right hand side of the equation against y
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T,
		double(T::*f)(double, double) const,
		double(T::*fx)(double, double) const,
		double(T::*fy)(double, double) const>
		double HigherOrderTaylorSeries(double x0, double y0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double x0_, y0_, y1, y1_;
			double h;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			if (x0 == x1)
			{
				return y0;
			}

			y1 = y0;
			h = (x1 - x0) / 10;

			x0_ = x0;
			y0_ = y0;

			do
			{
				y1_ = y1;
				while ((x1 - x0) * h > PRECISION)
				{
					y0 += h * (obj.*f)(x0, y0) + 0.5 * h * h * ((obj.*fx)(x0, y0)
							+ (obj.*fy)(x0, y0) * (obj.*f)(x0, y0));
					x0 += h;
				}
				y1 = y0;
				h /= 2;
				x0 = x0_;
				y0 = y0_;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

		/*!
		 * \brief Solve an ODE equation with Runge-Kutta method
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T, double(T::*f)(double, double) const>
		double RungeKutta(double x0, double y0, double x1,
				double gamma2, double tol, int maxIterations, const T& obj)
		{
			double gamma1, alpha, beta, h, x0_, y0_, y1, y1_, F_x_y_h;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			if (std::fabs(gamma2) < PRECISION)
			{
				gamma2 = 0.5;
			}

			alpha = 0.5 / gamma2;
			beta = 0.5 / gamma2;
			gamma1 = 1 - gamma2;

			y1 = y0;
			h = (x1 - x0) / 10;

			x0_ = x0;
			y0_ = y0;

			double alpha_mul_h, beta_mul_h;

			do
			{
				y1_ = y1;
				alpha_mul_h = alpha * h;
				beta_mul_h = beta * h;
				while ((x1 - x0) * h > PRECISION)
				{
					F_x_y_h = gamma1 * (obj.*f)(x0, y0)
							+ gamma2 * (obj.*f)(x0 + alpha_mul_h,
											y0 + beta_mul_h * (obj.*f)(x0, y0));
					y0 += h * F_x_y_h;
					x0 += h;
				}
				y1 = y0;
				h /= 2;
				x0 = x0_;
				y0 = y0_;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

		/*!
		 * \brief Solve an ODE equation with higher order Runge-Kutta method
		 *
		 * \tparam T The ODE object
		 * \tparam (T::*f) The right hand side of the equation
		 *
		 * \param x0 Initial value of x
		 * \param y0 Initial value of y
		 * \param x1 The ODE to be solved at
		 * \param tol Tolerance
		 * \param maxIterations Maximum number of iterations. This is to prevent
		 * non-converging problem and infinite loop. It is recommended to set 1000
		 * \param obj The equation
		 *
		 * \return The value y1, at x1
		 */
		template<class T, double(T::*f)(double, double) const>
		double RangeKuttaHigherOrder(double x0, double y0, double x1,
				double tol, int maxIterations, const T& obj)
		{
			double h, x0_, y0_, y1, y1_, F_x_y_h;
			int count = 0;

			if (maxIterations <= 0)
			{
				maxIterations = MAX_ITERATIONS;
			}

			y1 = y0;
			h = (x1 - x0) / 10;

			x0_ = x0;
			y0_ = y0;

			double v1, v2, v3, v4, h_2;

			do
			{
				y1_ = y1;
				h_2 = h / 2;
				while ((x1 - x0) * h > PRECISION)
				{
					v1 = (obj.*f)(x0, y0);
					v2 = (obj.*f)(x0 + h_2, y0 + h_2 * v1);
					v3 = (obj.*f)(x0 + h_2, y0 + h_2 * v2);
					v4 = (obj.*f)(x0 + h, y0 + h * v3);
					F_x_y_h = (v1 + 2 * v2 + 2 * v3 + v4) / 6;
					y0 += h * F_x_y_h;
					x0 += h;
				}
				y1 = y0;
				h /= 2;
				x0 = x0_;
				y0 = y0_;

			} while (std::fabs(y1_ - y1) > tol && ++count < maxIterations);

			return y1;
		}

	}
}

#endif /* ODE_H_ */
