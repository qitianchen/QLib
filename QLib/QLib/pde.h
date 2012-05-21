/*!
 * \file pde.h
 *
 * \brief A set of partial differential equation (PDE) solvers
 *
 * This file implements numeric solvers for two and three dimensional PDEs,
 * specified as follows, \n
 * Ut = a(t,x)*Uxx+b(t,x)*Ux+c(t,x)*U+d(t,x), and \n
 * Ut = a(t,x,y)*Uxx+b(t,x,y)*Uxy+c(t,x,y)*Uyy+
 * d(t,x,y)*Ux+e(t,x,y)*Uy+f(t,x,y)*U+g(t,x,y) \n
 *
 * They are suitable to solve Black-Scholes-Merton equations, for derivatives
 * of one underlining or two underlining assets.
 *
 *  Created on: May 6, 2012
 *      Author: Qitian Chen
 */

#ifndef PDE_H_
#define PDE_H_

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

#include "common.h"
#include "misc.h"

using namespace Eigen;

namespace QLib
{
	namespace PDE
	{
		#define DEBUG_FINAL_VALUES
		//#define DEBUG_STEP_VALUES

		/*!
		 * \brief Solve PDE in explicit method (forward euler)
		 * Ut = a(t,x)*Uxx+b(t,x)*Ux+c(t,x)*U+d(t,x)
		 *
		 * \tparam Equation PDE class
		 * \tparam (Equation::*a) Coefficient of Uxx
		 * \tparam (Equation::*b) Coefficient of Ux
		 * \tparam (Equation::*c) Coefficient of U
		 * \tparam (Equation::*d) Additional parameter
		 * \tparam InitialConditions Initial conditions
		 * \tparam BoundaryConditions Boundary conditions
		 * \tparam (BoundaryConditions::*lb) Values of U at lower boundary of x
		 * \tparam (BoundaryConditions::*ub) Values of U at upper boundary of x
		 *
		 * \param eqn The PDE problem object
		 * \param bc The boundary conditions object for x
		 * \param ic The initial conditions object for t
		 * \param t0 The initial timestamp
		 * \param t1 The terminal timestamp
		 * \param minX The lower bound of x
		 * \param maxX The upper bound of x
		 * \param Nt The number of time intervals
		 * \param Nx The number of space intervals
		 * \param x1 The terminal value of x
		 * \return U(t1,x1)
		 */
		template<class Equation,
		double (Equation::*a)(double, double) const,
		double (Equation::*b)(double, double) const,
		double (Equation::*c)(double, double) const,
		double (Equation::*d)(double, double) const,
		class InitialConditions,
		class BoundaryConditions,
		double (BoundaryConditions::*lb)(double, double) const,
		double (BoundaryConditions::*ub)(double, double) const>
		double Explicit(const Equation& eqn,
				const BoundaryConditions& bc,
				const InitialConditions& ic,
				double t0, double t1,
				double minX, double maxX,
				int Nt, int Nx, double x1)
		{
			VectorXd x = VectorXd::Zero(Nx+1);
			VectorXd v0 = VectorXd::Zero(Nx+1), v1 = v0;
			double k = (t1-t0)/Nt;
			double h = (maxX-minX)/Nx;
			int i, j;
			double a1, b1, c1, d1;

			// checking parameters
			assert(t0<t1);
			assert(minX<maxX);
			assert(minX<=x1);
			assert(x1<=maxX);

			for ( i = 0; i <= Nx; i++ ) // space
			{
				x[i] = i*h + minX;
				v0[i] = ic(x[i]);
			}

		#ifdef DEBUG_STEP_VALUES
			cout << "Lattice on space is " << endl << x << endl;
			cout << "Step found soln " << endl << v0 << endl;
		#endif

			// move forward time
			double t, A, B, C;
			double hq = h*h, h2 = 2*h, k_1 = 1/k;
			double tmp1, tmp2;
			t = t0;
			for ( j = 1; j <= Nt; j++ )
			{
				v1[0] = (bc.*lb)(t+k, minX);
				v1[Nx] = (bc.*ub)(t+k, maxX);
				for ( i = 1; i < Nx; i++ )
				{
					a1 = (eqn.*a)(t, x[i]);
					b1 = (eqn.*b)(t, x[i]);
					c1 = (eqn.*c)(t, x[i]);
					d1 = (eqn.*d)(t, x[i]);
					tmp1 = a1/hq; tmp2 = b1/h2;
					A = tmp1 - tmp2;
					B = -2*tmp1+c1+k_1;
					C = tmp1 + tmp2;
					v1[i] = k*(A*v0[i-1]+B*v0[i]+C*v0[i+1]+d1);
				}
				v0 = v1;
				t += k;
		#ifdef DEBUG_STEP_VALUES
				cout << "Step found soln " << endl << v0 << endl;
		#endif
			}

		#ifdef DEBUG_FINAL_VALUES
			cout << "Final explicit vector " << endl << v0 << endl;
		#endif

			i = round_double_to_int((x1-minX)/h);

			return v0[i];
		}

		/*!
		 * \brief Solve PDE in implicit method (backward euler)
		 * Ut = a(t,x)*Uxx+b(t,x)*Ux+c(t,x)*U+d(t,x)
		 *
		 * \tparam Equation PDE class
		 * \tparam (Equation::*a) Coefficient of Uxx
		 * \tparam (Equation::*b) Coefficient of Ux
		 * \tparam (Equation::*c) Coefficient of U
		 * \tparam (Equation::*d) Additional parameter
		 * \tparam InitialConditions Initial conditions
		 * \tparam BoundaryConditions Boundary conditions
		 * \tparam (BoundaryConditions::*lb) Values of U at lower boundary of x
		 * \tparam (BoundaryConditions::*ub) Values of U at upper boundary of x
		 *
		 * \param eqn The PDE problem object
		 * \param bc The boundary conditions object for x
		 * \param ic The initial conditions object for t
		 * \param t0 The initial timestamp
		 * \param t1 The terminal timestamp
		 * \param minX The lower bound of x
		 * \param maxX The upper bound of x
		 * \param Nt The number of time intervals
		 * \param Nx The number of space intervals
		 * \param x1 The terminal value of x
		 * \return U(t1,x1)
		 */
		template<class Equation,
		double (Equation::*a)(double, double) const,
		double (Equation::*b)(double, double) const,
		double (Equation::*c)(double, double) const,
		double (Equation::*d)(double, double) const,
		class InitialConditions,
		class BoundaryConditions,
		double (BoundaryConditions::*lb)(double, double) const,
		double (BoundaryConditions::*ub)(double, double) const>
		double Implicit(const Equation& eqn,
				const BoundaryConditions& bc,
				const InitialConditions& ic,
				double t0, double t1,
				double minX, double maxX,
				int Nt, int Nx, double x1)
		{
			VectorXd x = VectorXd::Zero(Nx+1);
			VectorXd v0 = VectorXd::Zero(Nx+1), v1 = v0;
			MatrixXd H = MatrixXd::Zero(Nx-1,Nx-1);
			MatrixXd L = MatrixXd::Zero(Nx-1,Nx-1);
			MatrixXd U = MatrixXd::Zero(Nx-1,Nx-1);
			VectorXd w = VectorXd::Zero(Nx-1);
			VectorXd v_ = VectorXd::Zero(Nx-1);
			double k = (t1-t0)/Nt;
			double h = (maxX-minX)/Nx;
			double a1, b1, c1, d1;
			int i, j;

			// checking parameters
			assert(t0<t1);
			assert(minX<maxX);
			assert(minX<=x1);
			assert(x1<=maxX);

			for ( i = 0; i <= Nx; i++ ) // space
			{
				x[i] = i*h + minX;
				v0[i] = ic(x[i]);
			}
			//cout << "Lattice on space is " << endl << x << endl;
			//cout << "Step found soln " << endl << v0 << endl;

			// move forward time
			double t, A, B, C, D;
			double hq = h*h, h2 = 2*h, k_1 = 1/k;
			double tmp1, tmp2;
			int ii;
			t = t0;
			for ( j = 1; j <= Nt; j++ )
			{
				t += k;
				v1[0] = (bc.*lb)(t, minX);
				v1[Nx] = (bc.*ub)(t, maxX);
				for ( i = 1; i < Nx; i++ )
				{
					a1 = (eqn.*a)(t, x[i]);
					b1 = (eqn.*b)(t, x[i]);
					c1 = (eqn.*c)(t, x[i]);
					d1 = (eqn.*d)(t, x[i]);
					tmp1 = a1/hq; tmp2 = b1/h2;
					A = tmp1 - tmp2;
					B = -2*tmp1+c1-k_1;
					C = tmp1 + tmp2;
					D = -d1 - v0[i]*k_1;
					ii = i-1;

					if ( i == 1 )
					{
						H(0,0) = B; H(0,1) = C;
						w[ii] = D - A*v1[0];
					}
					else if ( i == Nx-1 )
					{
						H(ii,ii-1) = A; H(ii,ii) = B;
						w[ii] = D - C*v1[Nx];
					}
					else
					{
						H(ii,ii-1) = A; H(ii,ii) = B; H(ii,ii+1) = C;
						w[ii] = D;
					}
				}
				//cout << "Implicit method step" << endl;
				//cout << "A is " << endl << H << endl;
				//cout << "b is " << endl << w << endl;
				Misc::TridiagonalLUFactorization(H, L, U);
				//cout << "LU decomposes A to " << endl;
				//cout << "L: " << endl << L << endl;
				//cout << "U: " << endl << U << endl;
				//cout << "L*U is " << endl << L*U << endl;
				v_ = U.inverse()*L.inverse()*w;
				//v_ = H.lu().solve(w);
				//cout << "A*x is " << endl << H*v_ << endl;
				v1.block(1,0,Nx-1,1) = v_;
				v0 = v1;
				//cout << "Step found soln " << endl << v0 << endl;
			}

		#ifdef DEBUG_FINAL_VALUES
			cout << "Final implicit vector " << endl << v0 << endl;
		#endif

			i = round_double_to_int((x1-minX)/h);

			return v0[i];
		}

		/*!
		 * \brief Solve PDE in Crank-Nicolson method
		 * Ut = a(t,x)*Uxx+b(t,x)*Ux+c(t,x)*U+d(t,x)
		 *
		 * \tparam Equation PDE class
		 * \tparam (Equation::*a) Coefficient of Uxx
		 * \tparam (Equation::*b) Coefficient of Ux
		 * \tparam (Equation::*c) Coefficient of U
		 * \tparam (Equation::*d) Additional parameter
		 * \tparam InitialConditions Initial conditions
		 * \tparam BoundaryConditions Boundary conditions
		 * \tparam (BoundaryConditions::*lb) Values of U at lower boundary of x
		 * \tparam (BoundaryConditions::*ub) Values of U at upper boundary of x
		 *
		 * \param eqn The PDE problem object
		 * \param bc The boundary conditions object for x
		 * \param ic The initial conditions object for t
		 * \param t0 The initial timestamp
		 * \param t1 The terminal timestamp
		 * \param minX The lower bound of x
		 * \param maxX The upper bound of x
		 * \param Nt The number of time intervals
		 * \param Nx The number of space intervals
		 * \param x1 The terminal value of x
		 * \return U(t1,x1)
		 */
		template<class Equation,
		double (Equation::*a)(double, double) const,
		double (Equation::*b)(double, double) const,
		double (Equation::*c)(double, double) const,
		double (Equation::*d)(double, double) const,
		class InitialConditions,
		class BoundaryConditions,
		double (BoundaryConditions::*lb)(double, double) const,
		double (BoundaryConditions::*ub)(double, double) const>
		double CrankNicolson(const Equation& eqn,
				const BoundaryConditions& bc,
				const InitialConditions& ic,
				double t0, double t1,
				double minX, double maxX,
				int Nt, int Nx, double x1)
		{
			VectorXd x = VectorXd::Zero(Nx+1);
			VectorXd v0 = VectorXd::Zero(Nx+1), v1 = v0;
			MatrixXd H = MatrixXd::Zero(Nx-1,Nx-1);
			MatrixXd L = MatrixXd::Zero(Nx-1,Nx-1);
			MatrixXd U = MatrixXd::Zero(Nx-1,Nx-1);
			VectorXd w = VectorXd::Zero(Nx-1);
			VectorXd v_ = VectorXd::Zero(Nx-1);
			double k = (t1-t0)/Nt;
			double h = (maxX-minX)/Nx;
			double a0, b0, c0, d0;
			double a1, b1, c1, d1;
			double t, t_, A, B, C, D;
			double A0, B0, C0;
			double hq = h*h, h2 = 2*h, k_1 = 1/k;
			double tmp1, tmp2;
			int ii;
			int i, j;

			// checking parameters
			assert(t0<t1);
			assert(minX<maxX);
			assert(minX<=x1);
			assert(x1<=maxX);

			for ( i = 0; i <= Nx; i++ ) // space
			{
				x[i] = i*h + minX;
				v0[i] = ic(x[i]);
			}

			t = t0;
			// move forward time
			for ( j = 1; j <= Nt; j++ )
			{
				t_ = t;
				t += k;
				v1[0] = (bc.*lb)(t, minX);
				v1[Nx] = (bc.*ub)(t, maxX);
				for ( i = 1; i < Nx; i++ )
				{
					a1 = (eqn.*a)(t, x[i]);
					b1 = (eqn.*b)(t, x[i]);
					c1 = (eqn.*c)(t, x[i]);
					d1 = (eqn.*d)(t, x[i]);
					tmp1 = a1/hq; tmp2 = b1/h2;
					A = tmp1 - tmp2;
					B = -2*tmp1+c1-k_1*2;
					C = tmp1 + tmp2;
					a0 = (eqn.*a)(t_, x[i]);
					b0 = (eqn.*b)(t_, x[i]);
					c0 = (eqn.*c)(t_, x[i]);
					d0 = (eqn.*d)(t_, x[i]);
					tmp1 = a0/hq; tmp2 = b0/h2;
					A0 = tmp1 - tmp2;
					B0 = -2*tmp1+c0+k_1*2;
					C0 = tmp1 + tmp2;
					D = -d0 - d1 - A0*v0[i-1] - B0*v0[i] - C0*v0[i+1];
					ii = i-1;

					if ( i == 1 )
					{
						H(0,0) = B; H(0,1) = C;
						w[ii] = D - A*v1[0];
					}
					else if ( i == Nx-1 )
					{
						H(ii,ii-1) = A; H(ii,ii) = B;
						w[ii] = D - C*v1[Nx];
					}
					else
					{
						H(ii,ii-1) = A; H(ii,ii) = B; H(ii,ii+1) = C;
						w[ii] = D;
					}
				}
				//cout << "Implicit method step" << endl;
				//cout << "A is " << endl << H << endl;
				//cout << "b is " << endl << w << endl;
				Misc::TridiagonalLUFactorization(H, L, U);
				//cout << "LU decomposes A to " << endl;
				//cout << "L: " << endl << L << endl;
				//cout << "U: " << endl << U << endl;
				//cout << "L*U is " << endl << L*U << endl;
				v_ = U.inverse()*L.inverse()*w;
				//v_ = H.lu().solve(w);
				//cout << "A*x is " << endl << H*v_ << endl;
				v1.block(1,0,Nx-1,1) = v_;
				v0 = v1;
			}

		#ifdef DEBUG_FINAL_VALUES
			cout << "Final Crank-Nicolson vector " << endl << v0 << endl;
		#endif

			i = round_double_to_int((x1-minX)/h);

			return v0[i];
		}

		/*!
		 * \brief Solve PDE in explicit method
		 * Ut = a(t,x,y)*Uxx+b(t,x,y)*Uxy+c(t,x,y)*Uyy+
		 * d(t,x,y)*Ux+e(t,x,y)*Uy+f(t,x,y)*U+g(t,x,y)
		 *
		 * \tparam Equation PDE class
		 * \tparam (Equation::*a) Coefficient of Uxx
		 * \tparam (Equation::*b) Coefficient of Uxy
		 * \tparam (Equation::*c) Coefficient of Uyy
		 * \tparam (Equation::*d) Coefficient of Ux
		 * \tparam (Equation::*e) Coefficient of Uy
		 * \tparam (Equation::*f) Coefficient of U
		 * \tparam (Equation::*g) Additional parameter
		 * \tparam InitialConditions Initial conditions
		 * \tparam BoundaryConditions Boundary conditions
		 * \tparam (BoundaryConditions::*lbx) Values of U at lower boundary of x
		 * \tparam (BoundaryConditions::*ubx) Values of U at upper boundary of x
		 * \tparam (BoundaryConditions::*lby) Values of U at lower boundary of y
		 * \tparam (BoundaryConditions::*ubt) Values of U at upper boundary of y
		 *
		 * \param eqn The PDE problem object
		 * \param bc The boundary conditions object for x and y
		 * \param ic The initial conditions object for t
		 * \param t0 The initial timestamp
		 * \param t1 The terminal timestamp
		 * \param minX The lower bound of x
		 * \param maxX The upper bound of x
		 * \param minY The lower bound of y
		 * \param maxY The upper bound of y
		 * \param Nt The number of time intervals
		 * \param Nx The number of space intervals for x
		 * \param Ny The number of space intervals for y
		 * \param x1 The terminal value of x
		 * \param y1 The terminal value of y
		 * \return U(t1,x1,y1)
		 */
		template<class Equation,
		double (Equation::*a)(double, double, double) const,
		double (Equation::*b)(double, double, double) const,
		double (Equation::*c)(double, double, double) const,
		double (Equation::*d)(double, double, double) const,
		double (Equation::*e)(double, double, double) const,
		double (Equation::*f)(double, double, double) const,
		double (Equation::*g)(double, double, double) const,
		class InitialConditions,
		class BoundaryConditions,
		double (BoundaryConditions::*lbx)(double, double, double) const,
		double (BoundaryConditions::*ubx)(double, double, double) const,
		double (BoundaryConditions::*lby)(double, double, double) const,
		double (BoundaryConditions::*uby)(double, double, double) const>
		double Explicit(const Equation& eqn,
				const BoundaryConditions& bc,
				const InitialConditions& ic,
				double t0, double t1,
				double minX, double maxX,
				double minY, double maxY,
				int Nt, int Nx, int Ny,
				double x1, double y1)
		{
			VectorXd x = VectorXd::Zero(Nx+1);
			VectorXd y = VectorXd::Zero(Ny+1);
			MatrixXd v0 = MatrixXd::Zero(Nx+1,Ny+1), v1 = v0;
			double k = (t1-t0)/Nt;
			double hx = (maxX-minX)/Nx;
			double hy = (maxY-minY)/Ny;
			int ix, iy, j;
			double a1, b1, c1, d1, e1, f1, g1;

			// checking parameters
			assert(t0<t1);
			assert(minX<maxX);
			assert(minY<maxY);
			assert(minX<=x1);
			assert(x1<=maxX);
			assert(minY<=y1);
			assert(y1<=maxY);

			for ( ix = 0; ix <= Nx; ix++ ) // space
			{
				x[ix] = ix*hx + minX;
			}

			for ( iy = 0; iy <= Ny; iy++ ) // space
			{
				y[iy] = iy*hy + minY;
			}

			for ( ix = 0; ix <= Nx; ix++ ) // space
			{
				for ( iy = 0; iy <= Ny; iy++ ) // space
				{
					v0(ix,iy) = ic(x[ix],y[iy]);
				}
			}
			//cout << "Lattice on space is " << endl << x << endl;
			//cout << "Step found soln " << endl << v0 << endl;

			// move forward time
			double t;
			double hxq = hx*hx, hx2 = 2*hx, k_1 = 1/k;
			double hyq = hy*hy, hy2 = 2*hy, hxy4 = 4*hx*hy;
			double tmp1, tmp2, tmp3, tmp4, tmp5;
			t = t0;
			for ( j = 1; j <= Nt; j++ )
			{
				t += k;
				// boundaries on x
				for ( iy = 0; iy <= Ny; iy++ )
				{
					v1(0,iy) = (bc.*lbx)(t, minX, y[iy]);
					v1(Nx,iy) = (bc.*ubx)(t, maxX, y[iy]);
				}
				// boundaries on y
				for ( ix = 0; ix <= Nx; ix++ )
				{
					v1(ix,0) = (bc.*lby)(t,x[ix], minY);
					v1(ix,Ny) = (bc.*uby)(t,x[ix], maxY);
				}
				t -= k;

				for ( ix = 1; ix < Nx; ix++ )
				{
					for ( iy = 1; iy < Ny; iy++ )
					{
						a1 = (eqn.*a)(t, x[ix], y[iy]);
						b1 = (eqn.*b)(t, x[ix], y[iy]);
						c1 = (eqn.*c)(t, x[ix], y[iy]);
						d1 = (eqn.*d)(t, x[ix], y[iy]);
						e1 = (eqn.*e)(t, x[ix], y[iy]);
						f1 = (eqn.*f)(t, x[ix], y[iy]);
						g1 = (eqn.*g)(t, x[ix], y[iy]);
						tmp1 = k*b1/(hxy4);
						tmp2 = a1/hxq;
						tmp3 = d1/hx2;
						tmp4 = c1/hyq;
						tmp5 = e1/hy2;
						v1(ix,iy) = tmp1*v0(ix-1,iy-1)+k*(tmp2-tmp3)*v0(ix-1,iy)-tmp1*v0(ix-1,iy+1)+
								k*(tmp4-tmp5)*v0(ix,iy-1)-k*(2*tmp2+2*tmp4-f1-k_1)*v0(ix,iy)+k*(tmp4+tmp5)*v0(ix,iy+1)-
								tmp1*v0(ix+1,iy-1)+k*(tmp2+tmp3)*v0(ix+1,iy)+tmp1*v0(ix+1,iy+1)+k*g1;
					}
				}
				v0 = v1;
				t += k;
				//cout << "Step found soln " << endl << v0 << endl;
			}

		#ifdef DEBUG_FINAL_VALUES
			cout << "Final explicit matrix " << endl << v0 << endl;
		#endif

			ix = round_double_to_int((x1-minX)/hx);
			iy = round_double_to_int((y1-minY)/hy);

			return v0(ix,iy);
		}

		/*!
		 * \brief Solve PDE in implicit method
		 * Ut = a(t,x,y)*Uxx+b(t,x,y)*Uxy+c(t,x,y)*Uyy+
		 * d(t,x,y)*Ux+e(t,x,y)*Uy+f(t,x,y)*U+g(t,x,y)
		 *
		 * \tparam Equation PDE class
		 * \tparam (Equation::*a) Coefficient of Uxx
		 * \tparam (Equation::*b) Coefficient of Uxy
		 * \tparam (Equation::*c) Coefficient of Uyy
		 * \tparam (Equation::*d) Coefficient of Ux
		 * \tparam (Equation::*e) Coefficient of Uy
		 * \tparam (Equation::*f) Coefficient of U
		 * \tparam (Equation::*g) Additional parameter
		 * \tparam InitialConditions Initial conditions
		 * \tparam BoundaryConditions Boundary conditions
		 * \tparam (BoundaryConditions::*lbx) Values of U at lower boundary of x
		 * \tparam (BoundaryConditions::*ubx) Values of U at upper boundary of x
		 * \tparam (BoundaryConditions::*lby) Values of U at lower boundary of y
		 * \tparam (BoundaryConditions::*ubt) Values of U at upper boundary of y
		 *
		 * \param eqn The PDE problem object
		 * \param bc The boundary conditions object for x and y
		 * \param ic The initial conditions object for t
		 * \param t0 The initial timestamp
		 * \param t1 The terminal timestamp
		 * \param minX The lower bound of x
		 * \param maxX The upper bound of x
		 * \param minY The lower bound of y
		 * \param maxY The upper bound of y
		 * \param Nt The number of time intervals
		 * \param Nx The number of space intervals for x
		 * \param Ny The number of space intervals for y
		 * \param x1 The terminal value of x
		 * \param y1 The terminal value of y
		 * \return U(t1,x1,y1)
		 */
		template<class Equation,
		double (Equation::*a)(double, double, double) const,
		double (Equation::*b)(double, double, double) const,
		double (Equation::*c)(double, double, double) const,
		double (Equation::*d)(double, double, double) const,
		double (Equation::*e)(double, double, double) const,
		double (Equation::*f)(double, double, double) const,
		double (Equation::*g)(double, double, double) const,
		class InitialConditions,
		class BoundaryConditions,
		double (BoundaryConditions::*lbx)(double, double, double) const,
		double (BoundaryConditions::*ubx)(double, double, double) const,
		double (BoundaryConditions::*lby)(double, double, double) const,
		double (BoundaryConditions::*uby)(double, double, double) const>
		double Implicit(const Equation& eqn,
				const BoundaryConditions& bc,
				const InitialConditions& ic,
				double t0, double t1,
				double minX, double maxX,
				double minY, double maxY,
				int Nt, int Nx, int Ny,
				double x1, double y1)
		{
			VectorXd x = VectorXd::Zero(Nx+1);
			VectorXd y = VectorXd::Zero(Ny+1);
			MatrixXd v0 = MatrixXd::Zero(Nx+1,Ny+1), v1 = v0;
			int N = (Nx-1)*(Ny-1);
			MatrixXd H = MatrixXd::Zero(N,N);
			MatrixXd L = H;
			MatrixXd U = H;
			VectorXd w = VectorXd::Zero(N);
			VectorXd v_ = w;
			MatrixXd Delta = MatrixXd::Zero(Ny-1,Ny-1), Sigma = Delta, Theta = Delta;
			VectorXd wi = VectorXd::Zero(Ny-1);
			double k = (t1-t0)/Nt;
			double hx = (maxX-minX)/Nx, hy = (maxY-minY)/Ny;
			double a1, b1, c1, d1, e1, f1, g1;
			int ix, iy, j;

			// checking parameters
			assert(t0<t1);
			assert(minX<maxX);
			assert(minY<maxY);
			assert(minX<=x1);
			assert(x1<=maxX);
			assert(minY<=y1);
			assert(y1<=maxY);

			for ( ix = 0; ix <= Nx; ix++ ) // space
			{
				x[ix] = ix*hx + minX;
			}

			for ( iy = 0; iy <= Ny; iy++ ) // space
			{
				y[iy] = iy*hy + minY;
			}

			for ( ix = 0; ix <= Nx; ix++ ) // space
			{
				for ( iy = 0; iy <= Ny; iy++ ) // space
				{
					v0(ix,iy) = ic(x[ix],y[iy]);
				}
			}

			//cout << "Lattice on space is " << endl << x << endl;
			//cout << "Step found soln " << endl << v0 << endl;

			// move forward time
			double t;
			double hxq = hx*hx, hx2 = 2*hx, k_1 = 1/k;
			double hyq = hy*hy, hy2 = 2*hy, hxy4 = 4*hx*hy;
			double tmp1, tmp2, tmp3, tmp4, tmp5;
			double A, B, C, D, E, F;
			int ix1, iy1, pos;
			t = t0;
			for ( j = 1; j <= Nt; j++ )
			{
				t += k;
				// boundaries on x
				for ( iy = 0; iy <= Ny; iy++ )
				{
					v1(0,iy) = (bc.*lbx)(t, minX, y[iy]);
					v1(Nx,iy) = (bc.*ubx)(t, maxX, y[iy]);
				}
				// boundaries on y
				for ( ix = 0; ix <= Nx; ix++ )
				{
					v1(ix,0) = (bc.*lby)(t,x[ix], minY);
					v1(ix,Ny) = (bc.*uby)(t,x[ix], maxY);
				}
				for ( ix = 1, ix1 = 0; ix < Nx; ix++, ix1++ )
				{
					for ( iy = 1, iy1 = 0; iy < Ny; iy++, iy1++ )
					{
						a1 = (eqn.*a)(t, x[ix], y[iy]);
						b1 = (eqn.*b)(t, x[ix], y[iy]);
						c1 = (eqn.*c)(t, x[ix], y[iy]);
						d1 = (eqn.*d)(t, x[ix], y[iy]);
						e1 = (eqn.*e)(t, x[ix], y[iy]);
						f1 = (eqn.*f)(t, x[ix], y[iy]);
						g1 = (eqn.*g)(t, x[ix], y[iy]);
						tmp1 = b1/(hxy4);
						tmp2 = a1/hxq;
						tmp3 = d1/hx2;
						tmp4 = c1/hyq;
						tmp5 = e1/hy2;
						A = tmp1;
						B = tmp2 + tmp3;
						C = tmp4 + tmp5;
						D = 2*(tmp2+tmp4)-f1+k_1;
						E = tmp4 - tmp5;
						F = tmp2 - tmp3;

						if ( ix == 1 )
						{
							wi[iy1] = -g1 - k_1*v0(ix,iy) - A*v1(ix1,iy1) - F*v1(ix1,iy) + A*v1(ix1,iy+1);
							if ( iy == 1 )
							{
								wi[iy1] += -E*v1(ix,iy1)+A*v1(ix+1,iy1);
							}
							else if ( iy == Ny-1 )
							{
								wi[iy1] -= C*v1(ix,Ny)+A*v1(ix+1,Ny);
							}
						}
						else if ( ix == Nx-1 )
						{
							wi[iy1] = -g1 - k_1*v0(ix,iy) - A*v1(ix1,iy1);
							if ( iy == 1 )
							{
								wi[iy1] -= E*v1(ix,iy1)-A*v1(Nx,iy1)+B*v1(Nx,iy)+A*v1(Nx,iy+1);
							}
							else if ( iy == Ny-1 )
							{
								wi[iy1] -= B*v1(Nx,iy)+A*v1(Nx,iy+1)-A*v1(ix-1,Ny)+C*v1(ix,Ny);
							}
							else
							{
								wi[iy1] -= B*v1(Nx,iy)+A*v1(Nx,iy+1);
							}
						}
						else
						{
							wi[iy1] = -g1 - k_1*v0(ix,iy);
							if ( iy == 1 )
							{
								wi[iy1] -= A*v1(ix1,0)+E*v1(ix,0)-A*v1(ix+1,0);
							}
							else if ( iy == Ny-1 )
							{
								wi[iy1] -= -A*v1(ix1,Ny)+C*v1(ix,Ny)+A*v1(ix+1,Ny);
							}
						}

						if ( iy == 1 )
						{
							Theta(iy1,iy1) = F; Theta(iy1,iy1+1) = -A;
							Delta(iy1,iy1) = -D; Delta(iy1,iy1+1) = C;
							Sigma(iy1,iy1) = B; Sigma(iy1,iy1+1) = A;
						}
						else if ( iy == Ny-1 )
						{
							Theta(iy1,iy1-1) = A; Theta(iy1,iy1) = F;
							Delta(iy1,iy1-1) = E; Delta(iy1,iy1) = -D;
							Sigma(iy1,iy1-1) = -A; Sigma(iy1,iy1) = B;
						}
						else
						{
							Theta(iy1,iy1-1) = A; Theta(iy1,iy1) = F; Theta(iy1,iy1+1) = -A;
							Delta(iy1,iy1-1) = E; Delta(iy1,iy1) = -D; Delta(iy1,iy1+1) = C;
							Sigma(iy1,iy1-1) = -A; Sigma(iy1,iy1) = B; Sigma(iy1,iy1+1) = A;
						}
					}

					pos = ix1*(Ny-1);
					w.block(pos,0,Ny-1,1) = wi;
					if ( ix == 1 )
					{
						H.block(pos,pos,Ny-1,Ny-1) = Delta;
						H.block(pos,pos+Ny-1,Ny-1,Ny-1) = Sigma;
					}
					else if ( ix == Nx-1 )
					{
						H.block(pos,pos-Ny+1,Ny-1,Ny-1) = Theta;
						H.block(pos,pos,Ny-1,Ny-1) = Delta;
					}
					else
					{
						H.block(pos,pos-Ny+1,Ny-1,Ny-1) = Theta;
						H.block(pos,pos,Ny-1,Ny-1) = Delta;
						H.block(pos,pos+Ny-1,Ny-1,Ny-1) = Sigma;
					}
				}
				v_ = Solver::LUFactorization(H, w);
				for ( ix = 1; ix < Nx; ix++ )
				{
		//			v1.block(ix,1,1,Ny-1) = v_.block(ix*(Ny-1),0,Ny-1,1).transpose();
					for ( iy = 1; iy < Ny; iy++ )
					{
						v1(ix,iy) = v_[(ix-1)*(Ny-1)+iy-1];
					}
				}
				v0 = v1;
			}

		#ifdef DEBUG_FINAL_VALUES
			cout << "Final implicit matrix " << endl << v0 << endl;
		#endif

			ix = round_double_to_int((x1-minX)/hx);
			iy = round_double_to_int((y1-minY)/hy);

			return v0(ix,iy);
		}


		/*!
		 * \brief Solve PDE in operator splitting method
		 * Ut = a(t,x,y)*Uxx+b(t,x,y)*Uxy+c(t,x,y)*Uyy+
		 * d(t,x,y)*Ux+e(t,x,y)*Uy+f(t,x,y)*U+g(t,x,y)
		 *
		 *
		 * \tparam Equation PDE class
		 * \tparam (Equation::*a) Coefficient of Uxx
		 * \tparam (Equation::*b) Coefficient of Uxy
		 * \tparam (Equation::*c) Coefficient of Uyy
		 * \tparam (Equation::*d) Coefficient of Ux
		 * \tparam (Equation::*e) Coefficient of Uy
		 * \tparam (Equation::*f) Coefficient of U
		 * \tparam (Equation::*g) Additional parameter
		 * \tparam InitialConditions Initial conditions
		 * \tparam BoundaryConditions Boundary conditions
		 * \tparam (BoundaryConditions::*lbx) Values of U at lower boundary of x
		 * \tparam (BoundaryConditions::*ubx) Values of U at upper boundary of x
		 * \tparam (BoundaryConditions::*lby) Values of U at lower boundary of y
		 * \tparam (BoundaryConditions::*ubt) Values of U at upper boundary of y
		 *
		 * \param eqn The PDE problem object
		 * \param bc The boundary conditions object for x and y
		 * \param ic The initial conditions object for t
		 * \param t0 The initial timestamp
		 * \param t1 The terminal timestamp
		 * \param minX The lower bound of x
		 * \param maxX The upper bound of x
		 * \param minY The lower bound of y
		 * \param maxY The upper bound of y
		 * \param Nt The number of time intervals
		 * \param Nx The number of space intervals for x
		 * \param Ny The number of space intervals for y
		 * \param x1 The terminal value of x
		 * \param y1 The terminal value of y
		 * \return U(t1,x1,y1)
		 */
		template<class Equation,
		double (Equation::*a)(double, double, double) const,
		double (Equation::*b)(double, double, double) const,
		double (Equation::*c)(double, double, double) const,
		double (Equation::*d)(double, double, double) const,
		double (Equation::*e)(double, double, double) const,
		double (Equation::*f)(double, double, double) const,
		double (Equation::*g)(double, double, double) const,
		class InitialConditions,
		class BoundaryConditions,
		double (BoundaryConditions::*lbx)(double, double, double) const,
		double (BoundaryConditions::*ubx)(double, double, double) const,
		double (BoundaryConditions::*lby)(double, double, double) const,
		double (BoundaryConditions::*uby)(double, double, double) const>
		double OperatorSplitting(const Equation& eqn,
				const BoundaryConditions& bc,
				const InitialConditions& ic,
				double t0, double t1,
				double minX, double maxX,
				double minY, double maxY,
				int Nt, int Nx, int Ny,
				double x1, double y1)
		{
			VectorXd x = VectorXd::Zero(Nx+1);
			VectorXd y = VectorXd::Zero(Ny+1);
			MatrixXd v0 = MatrixXd::Zero(Nx+1,Ny+1), v1 = v0, v2 = v0;
			MatrixXd Hx = MatrixXd::Zero(Nx-1,Nx-1), Hy = MatrixXd::Zero(Ny-1,Ny-1);
			MatrixXd Lx = Hx, Ly = Hy;
			MatrixXd Ux = Hx, Uy = Hy;
			VectorXd wx = VectorXd::Zero(Nx-1), wy = VectorXd::Zero(Ny-1);
			VectorXd v1_ = VectorXd::Zero(Nx-1), v2_ = VectorXd::Zero(Ny-1);
			double k = (t1-t0)/Nt;
			double hx = (maxX-minX)/Nx, hy = (maxY-minY)/Ny;
			double b0;
			double a1, b1, c1, d1, e1, f1, g1;
			int ix, iy, j;

			// checking parameters
			assert(t0<t1);
			assert(minX<maxX);
			assert(minY<maxY);
			assert(minX<=x1);
			assert(x1<=maxX);
			assert(minY<=y1);
			assert(y1<=maxY);

			for ( ix = 0; ix <= Nx; ix++ ) // space
			{
				x[ix] = ix*hx + minX;
			}

			for ( iy = 0; iy <= Ny; iy++ ) // space
			{
				y[iy] = iy*hy + minY;
			}

			for ( ix = 0; ix <= Nx; ix++ ) // space
			{
				for ( iy = 0; iy <= Ny; iy++ ) // space
				{
					v0(ix,iy) = ic(x[ix],y[iy]);
				}
			}

			//cout << "Lattice on space is " << endl << x << endl;
			//cout << "Step found soln " << endl << v0 << endl;

			// move forward time
			double t, t_;
			double hxq = hx*hx, hx2 = 2*hx, k_1 = 1/k;
			double hyq = hy*hy, hy2 = 2*hy, hxy4 = 4*hx*hy;
			double tmp1, tmp2;
			double A, B, C, F;
			int ix1, iy1;
			t = t0; t_ = t0;
			for ( j = 1; j <= Nt; j++, t_ += k )
			{
				t += k;
				// boundaries on x
				for ( iy = 0; iy <= Ny; iy++ )
				{
					v1(0,iy) = (bc.*lbx)(t, minX, y[iy]);
					v1(Nx,iy) = (bc.*ubx)(t, maxX, y[iy]);
					v2(0,iy) = (bc.*lbx)(t, minX, y[iy]);
					v2(Nx,iy) = (bc.*ubx)(t, maxX, y[iy]);
				}
				// boundaries on y
				for ( ix = 0; ix <= Nx; ix++ )
				{
					v1(ix,0) = (bc.*lby)(t,x[ix], minY);
					v1(ix,Ny) = (bc.*uby)(t,x[ix], maxY);
					v2(ix,0) = (bc.*lby)(t,x[ix], minY);
					v2(ix,Ny) = (bc.*uby)(t,x[ix], maxY);
				}
				// step 1
				for ( iy = 1, iy1 = 0; iy < Ny; iy++, iy1++ )
				{
					for ( ix = 1, ix1 = 0; ix < Nx; ix++, ix1++ )
					{
						b0 = (eqn.*b)(t_, x[ix], y[iy]);
						a1 = (eqn.*a)(t, x[ix], y[iy]);
						b1 = (eqn.*b)(t, x[ix], y[iy]);
						d1 = (eqn.*d)(t, x[ix], y[iy]);
						f1 = (eqn.*f)(t, x[ix], y[iy]);
						g1 = (eqn.*g)(t, x[ix], y[iy]);
						tmp1 = a1/hxq;
						tmp2 = d1/hx2;
						A = tmp1-tmp2;
						B = -(2*tmp1+k_1-f1);
						C = tmp1+tmp2;
						F= -0.5*b0*(v0(ix+1,iy+1)-v0(ix+1,iy-1)-
								v0(ix-1,iy+1)+v0(ix-1,iy-1))/hxy4-
								k_1*v0(ix,iy)-g1;
						if ( ix == 1 )
						{
							Hx(ix1,ix1) = B; Hx(ix1,ix) = C;
							wx[ix1] = F - A*v1(ix1,iy);
						}
						else if ( ix == Nx-1 )
						{
							Hx(ix1,ix1-1) = A; Hx(ix1,ix1) = B;
							wx[ix1] = F - C*v1(Nx,iy);
						}
						else
						{
							Hx(ix1,ix1-1) = A; Hx(ix1,ix1) = B; Hx(ix1,ix) = C;
							wx[ix1] = F;
						}
					}
					Misc::TridiagonalLUFactorization(Hx, Lx, Ux);
					v1_ = Ux.inverse()*Lx.inverse()*wx;
					v1.block(1,iy,Nx-1,1) = v1_;
				}
				// step 2
				for ( ix = 1, ix1 = 0; ix < Nx; ix++, ix1++ )
				{
					for ( iy = 1, iy1 = 0; iy < Ny; iy++, iy1++ )
					{
						b0 = (eqn.*b)(t_, x[ix], y[iy]);
						b1 = (eqn.*b)(t, x[ix], y[iy]);
						c1 = (eqn.*c)(t, x[ix], y[iy]);
						e1 = (eqn.*e)(t, x[ix], y[iy]);
						tmp1 = c1/hyq;
						tmp2 = e1/hy2;
						A = tmp1-tmp2;
						B = -(2*tmp1+k_1);
						C = tmp1+tmp2;
						F= -0.5*b0*(v1(ix+1,iy+1)-v1(ix+1,iy-1)-
								v1(ix-1,iy+1)+v1(ix-1,iy-1))/hxy4-
								k_1*v1(ix,iy);
						if ( iy == 1 )
						{
							Hy(iy1,iy1) = B; Hy(iy1,iy) = C;
							wy[iy1] = F - A*v2(ix,iy1);
						}
						else if ( iy == Ny-1 )
						{
							Hy(iy1,iy1-1) = A; Hy(iy1,iy1) = B;
							wy[iy1] = F - C*v2(ix,Ny);
						}
						else
						{
							Hy(iy1,iy1-1) = A; Hy(iy1,iy1) = B; Hy(iy1,iy) = C;
							wy[iy1] = F;
						}
					}
					Misc::TridiagonalLUFactorization(Hy, Ly, Uy);
					v2_ = Uy.inverse()*Ly.inverse()*wy;
					v2.block(ix,1,1,Ny-1) = v2_.transpose();
				}
				v0 = v2;
			}

		#ifdef DEBUG_FINAL_VALUES
			cout << "Final operator splitting matrix " << endl << v0 << endl;
		#endif

			ix = round_double_to_int((x1-minX)/hx);
			iy = round_double_to_int((y1-minY)/hy);

			return v0(ix,iy);
		}

	}
}

#endif /* PDE_H_ */
