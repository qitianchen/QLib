/*
 * BlackScholesMertonEquation.cpp
 *
 *  Created on: May 6, 2012
 *      Author: Administrator
 */

#include "BlackScholesMertonEquation.h"

BlackScholesMertonEquation::BlackScholesMertonEquation(double r, double vol)
:r(r), vol(vol)
{
	// TODO Auto-generated constructor stub
	a0 = 0.5*vol*vol;
	b0 = r;
	c0 = -r;
}

BlackScholesMertonEquation::~BlackScholesMertonEquation()
{
	// TODO Auto-generated destructor stub
}

double BlackScholesMertonEquation::a(double t, double S) const
{
	return a0*S*S;
}

double BlackScholesMertonEquation::b(double t, double S) const
{
	return b0*S;
}

double BlackScholesMertonEquation::c(double t, double S) const
{
	return c0;
}

double BlackScholesMertonEquation::d(double t, double S) const
{
	return 0;
}
