/*
 * BlackScholesMertonEquation.h
 *
 *  Created on: May 6, 2012
 *      Author: Administrator
 */

#ifndef BLACKSCHOLESMERTONEQUATION_H_
#define BLACKSCHOLESMERTONEQUATION_H_

class BlackScholesMertonEquation
{
public:
	BlackScholesMertonEquation(double r, double vol);
	virtual ~BlackScholesMertonEquation();
	virtual double a(double t, double S) const;
	virtual double b(double t, double S) const;
	virtual double c(double t, double S) const;
	virtual double d(double t, double S) const;

private:
	double r, vol;
	double a0, b0, c0;
};

#endif /* BLACKSCHOLESMERTONEQUATION_H_ */
