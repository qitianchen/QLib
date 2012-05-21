//============================================================================
// Name        : QLib.cpp
// Author      : Chen Qitian
// Version     :
// Copyright   :
// Description : QLib Examples
//============================================================================

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "QLib/solver.h"
#include "QLib/solver_linear_system.h"
#include "QLib/ode.h"
#include "QLib/pde.h"
#include "QLib/misc.h"

#include "BlackScholesMertonEquation.h"
#include "PayOffCall.h"
#include "BoundaryCall.h"
#include "BarrierPayoff.h"
#include "BarrierBoundaries.h"
#include "TwoFactorsBlackScholesMertonEquation.h"
#include "ExchangePutPayoff.h"
#include "ExchangePutBoundaries.h"

using namespace std;
using namespace Eigen;
using namespace QLib;

/*!
 * \brief A polynomial equation
 * x^6-x-1=0
 */
class NonLinearEquation
{
public:
	inline double f(double x) const
	{
		return pow(x, 6) - x - 1;
	}
	inline double operator ()(double x) const
	{
		return f(x);
	}
	inline double fx(double x) const
	{
		return 6*pow(x, 5) - 1;
	}
	inline double g(double x) const
	{
		return pow(x+1, 1.0/6);
	}
};

/*!
 * \brief Test cases of non-linear equations
 */
void test_non_linear_equation()
{
	NonLinearEquation eqn;
	double x0 = 0, x1 = 1, x2 = 2;

	// Bisection solver, no operator overloading
	x0 = Solver::Bisection<NonLinearEquation, &NonLinearEquation::f>
		(0, x1, x2, 1e-5, 1000, eqn);
	cout<<"Root from Bisection 1: "<<x0<<endl;

	// Bisection solver, operator overloading
	x0 = Solver::Bisection<NonLinearEquation>(0, x1, x2, 1e-5, 1000, eqn);
	cout<<"Root from Bisection 2: "<<x0<<endl;

	// Newton-Raphson method
	x0 = Solver::NewtonRaphson<NonLinearEquation,
			&NonLinearEquation::f,
			&NonLinearEquation::fx>(0, x0, 1e-5, 1000, eqn);
	cout<<"Root from Newton-Raphson: "<<x0<<endl;

	// Secant method, no operator overloading
	x0 = Solver::Secant<NonLinearEquation,
			&NonLinearEquation::f>(0, x1, x2, 1e-5, 1000, eqn);
	cout<<"Root from Secant 1: "<<x0<<endl;

	// Secant method, operator overloading
	x0 = Solver::Secant<NonLinearEquation>(0, x1, x2, 1e-5, 1000, eqn);
	cout<<"Root from Secant 2: "<<x0<<endl;

	// Fixed point iteration
	x0 = Solver::FixedPointIteration<NonLinearEquation,
			&NonLinearEquation::f,
			&NonLinearEquation::g>(0, x0, 1e-5, 1000, eqn);
	cout<<"Root from Fixed Point Iteration: "<<x0<<endl;
}

/*!
 * \brief Test Cases for system of linear equations
 */
void test_linear_system()
{
	Matrix4d A;
	A << 4,3,2,1,
		3,4,3,2,
		2,3,4,3,
		1,2,3,4;
	Vector4d b;
	b << 1,1,-1,-1 ;
	Vector4d x;
	x << 0,0,0,0;

	cout << "The problem is A = " << endl << A << endl;
	cout << "b = " << endl << b << endl;

	// inverse matrix on the left
	x = Solver::DirectInversion(A, b);
	cout << "The solution from direct inversion is " <<endl << x << endl;
	cout << "Verify: A*x = " << endl << A*x << endl;

	// Gaussian elimination
	x = Solver::GaussianElimination(A, b);
	cout << "The solution from Gaussian Elimination is " <<endl << x << endl;
	cout << "Verify: A*x = " << endl << A*x << endl;

	// Jacobi Iterative method
	VectorXd x0 = VectorXd::Zero(4);
	x = Solver::JacobiIterative(A, b, x0, 1e-5, 0);
	cout << "The solution from Jacobi Iterative is " <<endl << x << endl;
	cout << "Verify: A*x = " << endl << A*x << endl;

	// Gauss-Seibal method
	x = Solver::GaussSeidal(A, b, x0, 1e-5, 0);
	cout << "The solution from Gauss-Seidal is " <<endl << x << endl;
	cout << "Verify: A*x = " << endl << A*x << endl;

	// LU facotorization
	x = Solver::LUFactorization(A, b);
	cout << "The solution from LU Factorization is " <<endl << x << endl;
	cout << "Verify: A*x = " << endl << A*x << endl;
}

/*!
 * \brief ODE
 * y'(x) = -2*x*y+2*x
 */
class OrdinaryDifferentialEquation
{
public:
	inline double f(double x, double y) const
	{
		return -2*x*y+2*x;
	}
	inline double fx(double x, double y) const
	{
		return -2*y+2;
	}
	inline double fy(double x, double y) const
	{
		return -2*x;
	}
	inline double solver1(double x, double y, double h) const
	{
		double tmp = 2*h*x;
		return (y+tmp)/(1+tmp);
	}
};
/*!
 * \brief Test cases for Ordinary Differential Equations
 */
void test_ode()
{
	OrdinaryDifferentialEquation testode;
	double x0 = 0, y0 = 2, x1 = 1.5;

	double y1 = ODE::Difference<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::f>(x0, y0, x1, 1e-5, 0, testode);
	cout<<"Solution from Forward Difference: "<<y1<<endl;

	y1 = ODE::RichardsonExtrapolation<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::f>(x0, y0, x1, 1e-5, 1000, testode);
	cout<<"Solution from Richardson Extrapolation: "<<y1<<endl;

	y1 = ODE::Implicit<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::solver1>(x0, y0, x1, 1e-5, 1000, testode);
	cout<<"Solution from Implicit method: "<<y1<<endl;

	y1 = ODE::PredictorCorrector<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::f>(x0, y0, x1, 1e-5, 1000, testode);
	cout<<"Solution from Predictor Corrector method: "<<y1<<endl;

	y1 = ODE::Trapezium<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::f>(x0, y0, x1, 1e-5, 1000, testode);
	cout<<"Solution from Trapezium method: "<<y1<<endl;

	y1 = ODE::HigherOrderTaylorSeries<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::f,
			&OrdinaryDifferentialEquation::fx,
			&OrdinaryDifferentialEquation::fy>(x0, y0, x1, 1e-5, 1000, testode);
	cout<<"Solution from Higher Order Taylor Series method: "<<y1<<endl;

	y1 = ODE::RungeKutta<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::f>(x0, y0, x1, 1, 1e-5, 1000, testode);
	cout<<"Solution from Runge-Kutta method: "<<y1<<endl;

	y1 = ODE::RangeKuttaHigherOrder<OrdinaryDifferentialEquation,
			&OrdinaryDifferentialEquation::f>(x0, y0, x1, 1e-5, 1000, testode);
	cout<<"Solution from Higher Order Runge-Kutta method: "<<y1<<endl;
}


/*!
 * \brief A generic 2-D PDE equation
 * Ut = a(t,x)*Uxx+b(t,x)*Ux+c(t,x)*U+d(t,x)
 */
class PartialDifferentialEquation
{
public:
	PartialDifferentialEquation(double a, double b, double c, double d)
		:a0(a), b0(b), c0(c), d0(d){};
	virtual ~PartialDifferentialEquation(){};
	virtual double a(double t, double x) const
	{
		return a0;
	}
	virtual double b(double t, double x) const
	{
		return b0;
	}
	virtual double c(double t, double x) const
	{
		return c0;
	}
	virtual double d(double t, double x) const
	{
		return d0;
	}

private:
	double a0, b0, c0, d0;
};

/*!
 * \brief Boundary conditions
 * Lower bound v(t,lb) = exp(t)
 * Upper bound v(t,ub) = exp(1+t)
 */
class TestBoundaryConditions
{
public:
	TestBoundaryConditions(double lb, double ub):lb0(lb),ub0(ub){};
	double lb(double t, double x) const
	{
		return std::exp(t);
	};
	double ub(double t, double x) const
	{
		return std::exp(t+1);
	};

private:
	double lb0, ub0;
};

/*!
 * \brief Initial conditions
 * v(t0,x) = exp(x)
 */
class TestInitialConditions
{
public:
	double operator ()(double x) const
	{
		return std::exp(x);
	};
};

/*!
 * \brief Test cases for two dimensional partial differential equations
 */
void test_pde()
{
	PartialDifferentialEquation mypde(1, -1, 1, 0);
	TestBoundaryConditions mybc(0, 1);
	TestInitialConditions myic;

	double v = PDE::Explicit<PartialDifferentialEquation,
			&PartialDifferentialEquation::a,
			&PartialDifferentialEquation::b,
			&PartialDifferentialEquation::c,
			&PartialDifferentialEquation::d,
			TestInitialConditions,
			TestBoundaryConditions,
			&TestBoundaryConditions::lb,
			&TestBoundaryConditions::ub>(mypde, mybc, myic,
					0, 1, 0, 1,
					4, 4, // this is not stable, only for illustration
					0.5);
	cout << "Solution from Explicit method is " << v << endl;

	v = PDE::Implicit<PartialDifferentialEquation,
			&PartialDifferentialEquation::a,
			&PartialDifferentialEquation::b,
			&PartialDifferentialEquation::c,
			&PartialDifferentialEquation::d,
			TestInitialConditions,
			TestBoundaryConditions,
			&TestBoundaryConditions::lb,
			&TestBoundaryConditions::ub>(mypde, mybc, myic,
					0, 1, 0, 1,
					4, 4,
					0.5);
	cout << "Solution from Implicit method is " << v << endl;

	v = PDE::CrankNicolson<PartialDifferentialEquation,
			&PartialDifferentialEquation::a,
			&PartialDifferentialEquation::b,
			&PartialDifferentialEquation::c,
			&PartialDifferentialEquation::d,
			TestInitialConditions,
			TestBoundaryConditions,
			&TestBoundaryConditions::lb,
			&TestBoundaryConditions::ub>(mypde, mybc, myic,
					0, 1, 0, 1,
					4, 4,
					0.5);
	cout << "Solution from Crank-Nicolson method is " << v << endl;
}

/*!
 * \brief Test cases for one risk factor Black-Scholes-Merton equations
 */
void test_bsm_pde()
{
	// European call option
	// interest rate 6%, volatility 20%
	// strike price $50
	BlackScholesMertonEquation bsm(0.06, 0.2);
	PayOffCall call(50);
	// strike price and interest rate for boundary conditions
	BoundaryCall boundaries(50, 0.06);

	// explicit method
	double v = PDE::Explicit<BlackScholesMertonEquation,
			&BlackScholesMertonEquation::a, // coefficient of Vss
			&BlackScholesMertonEquation::b, // coefficient of Vs
			&BlackScholesMertonEquation::c, // coefficient of V
			&BlackScholesMertonEquation::d, // dummy
			PayOffCall,
			BoundaryCall,
			&BoundaryCall::lb,
			&BoundaryCall::ub>(bsm, boundaries, call,
					0, 0.4, 0, 100, // range
					10, 10, // number of intervals in time and space
					50); // S0
	cout << "Call option price from Explicit method is " << v << endl;

	// implicit method
	v = PDE::Implicit<BlackScholesMertonEquation,
			&BlackScholesMertonEquation::a,
			&BlackScholesMertonEquation::b,
			&BlackScholesMertonEquation::c,
			&BlackScholesMertonEquation::d,
			PayOffCall,
			BoundaryCall,
			&BoundaryCall::lb,
			&BoundaryCall::ub>(bsm, boundaries, call,
					0, 0.4, 0, 100,
					10, 10,
					50);
	cout << "Call option price from Implicit method is " << v << endl;

	// Crank-Nicolson method
	v = PDE::CrankNicolson<BlackScholesMertonEquation,
			&BlackScholesMertonEquation::a,
			&BlackScholesMertonEquation::b,
			&BlackScholesMertonEquation::c,
			&BlackScholesMertonEquation::d,
			PayOffCall,
			BoundaryCall,
			&BoundaryCall::lb,
			&BoundaryCall::ub>(bsm, boundaries, call,
					0, 0.4, 0, 100,
					10, 10,
					50);
	cout << "Call option price from Crank-Nicolson method is " << v << endl;

	// Barrier option
	// interest rate 10%, volatility 20%
	// strike price $100, lower barrier $75, upper barrier $130
	BlackScholesMertonEquation bsm2(0.1, 0.2);
	BarrierPayoff barrier(100, 75, 130);
	BarrierBoundaries barriers(75, 130);
	v = PDE::CrankNicolson<BlackScholesMertonEquation,
			&BlackScholesMertonEquation::a,
			&BlackScholesMertonEquation::b,
			&BlackScholesMertonEquation::c,
			&BlackScholesMertonEquation::d,
			BarrierPayoff,
			BarrierBoundaries,
			&BarrierBoundaries::lb,
			&BarrierBoundaries::ub>(bsm2, barriers, barrier,
					0, 1, 75, 130, // range
					10, 10, // number of intervals in time and space
					90); //S0
	cout << "Barrier option price from Crank-Nicolson method is " << v << endl;

}

/*!
 * \brief Test cases for two risk factors Black-Scholes-Merton equations
 */
void test_bsm_pde_2_factors()
{
	// exchange put option
	// interest rate 3%, volatility of asset 1 20%, asset 2 10%, correlation 0.05
	TwoFactorsBlackScholesMertonEquation bsm2f(0.03, 0.2, 0.1, 0.05);
	ExchangePutPayoff exchangePutOption;
	ExchangePutBoundaries exchangeBoundaries(100, 100);

	double v = PDE::Explicit<TwoFactorsBlackScholesMertonEquation,
			&TwoFactorsBlackScholesMertonEquation::a, // coefficient of Vs1s1
			&TwoFactorsBlackScholesMertonEquation::b, // coefficient of Vs1s2
			&TwoFactorsBlackScholesMertonEquation::c, // coefficient of Vs2s2
			&TwoFactorsBlackScholesMertonEquation::d, // coefficient of Vs1
			&TwoFactorsBlackScholesMertonEquation::e, // coefficient of Vs2
			&TwoFactorsBlackScholesMertonEquation::f, // coefficient of V
			&TwoFactorsBlackScholesMertonEquation::g, // dummy
			ExchangePutPayoff,
			ExchangePutBoundaries,
			&ExchangePutBoundaries::lbx,
			&ExchangePutBoundaries::ubx,
			&ExchangePutBoundaries::lby,
			&ExchangePutBoundaries::uby>(bsm2f, exchangeBoundaries, exchangePutOption,
					0, 0.8, 0, 100, 0, 100, // range
					100, 10, 10, // number of intervals
					50, 50); // (S1_0, S2_0)
	cout.width(10);
	cout << "Exchange option price from Explicit method is " << v << endl;

	v = PDE::Implicit<TwoFactorsBlackScholesMertonEquation,
			&TwoFactorsBlackScholesMertonEquation::a,
			&TwoFactorsBlackScholesMertonEquation::b,
			&TwoFactorsBlackScholesMertonEquation::c,
			&TwoFactorsBlackScholesMertonEquation::d,
			&TwoFactorsBlackScholesMertonEquation::e,
			&TwoFactorsBlackScholesMertonEquation::f,
			&TwoFactorsBlackScholesMertonEquation::g,
			ExchangePutPayoff,
			ExchangePutBoundaries,
			&ExchangePutBoundaries::lbx,
			&ExchangePutBoundaries::ubx,
			&ExchangePutBoundaries::lby,
			&ExchangePutBoundaries::uby>(bsm2f, exchangeBoundaries, exchangePutOption,
					0, 0.8, 0, 100, 0, 100,
					100, 10, 10,
					50, 50);
	cout.width(10);
	cout << "Exchange option price from Implicit method is " << v << endl;

	v = PDE::OperatorSplitting<TwoFactorsBlackScholesMertonEquation,
			&TwoFactorsBlackScholesMertonEquation::a,
			&TwoFactorsBlackScholesMertonEquation::b,
			&TwoFactorsBlackScholesMertonEquation::c,
			&TwoFactorsBlackScholesMertonEquation::d,
			&TwoFactorsBlackScholesMertonEquation::e,
			&TwoFactorsBlackScholesMertonEquation::f,
			&TwoFactorsBlackScholesMertonEquation::g,
			ExchangePutPayoff,
			ExchangePutBoundaries,
			&ExchangePutBoundaries::lbx,
			&ExchangePutBoundaries::ubx,
			&ExchangePutBoundaries::lby,
			&ExchangePutBoundaries::uby>(bsm2f, exchangeBoundaries, exchangePutOption,
					0, 0.8, 0, 100, 0, 100,
					100, 10, 10,
					50, 50);
	cout.width(10);
	cout << "Exchange option price from Operator Splitting method is " << v << endl;
}

int main() {
	cout << "Hello QLib" << endl; // prints

	cout << "=========================================" << endl << endl;
	cout << "Test non-linear solver" << endl;
	test_non_linear_equation();
	cout << "=========================================" << endl << endl;
	cout << "Test system of linear equations " << endl;
	test_linear_system();
	cout << "=========================================" << endl << endl;
	cout << "Test Ordinary Differential Equation " << endl;
	test_ode();
	cout << "=========================================" << endl << endl;
	cout << "Test Partial Differential Equation " << endl;
	test_pde();
	cout << "=========================================" << endl << endl;
	cout << "Test One Factor Black-Scholes-Merton Equation " << endl;
	test_bsm_pde();
	cout << "=========================================" << endl << endl;
	cout << "Test Two Factors Black-Scholes-Merton Equation " << endl;
	test_bsm_pde_2_factors();
	cout << "=========================================" << endl << endl;

	return 0;
}
