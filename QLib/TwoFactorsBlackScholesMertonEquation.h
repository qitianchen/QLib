/*
 * TwoFactorsBlackScholesMertonEquation.h
 *
 *  Created on: May 14, 2012
 *      Author: Administrator
 */

#ifndef TWOFACTORSBLACKSCHOLESMERTONEQUATION_H_
#define TWOFACTORSBLACKSCHOLESMERTONEQUATION_H_

class TwoFactorsBlackScholesMertonEquation
{
public:
	TwoFactorsBlackScholesMertonEquation(double r, double vol1, double vol2, double correl);
	virtual ~TwoFactorsBlackScholesMertonEquation();
	inline double a(double t, double S1, double S2) const
	{
		return tmp1*S1*S1;
	}
	inline double b(double t, double S1, double S2) const
	{
		return tmp3*S1*S2;
	}
	inline double c(double t, double S1, double S2) const
	{
		return tmp2*S2*S2;
	}
	inline double d(double t, double S1, double S2) const
	{
		return r*S1;
	}
	inline double e(double t, double S1, double S2) const
	{
		return r*S2;
	}
	inline double f(double t, double S1, double S2) const
	{
		return -r;
	}
	inline double g(double t, double S1, double S2) const
	{
		return 0;
	}

private:
	double r;
	double vol1, vol2;
	double correl;
	double tmp1, tmp2, tmp3;
};

#endif /* TWOFACTORSBLACKSCHOLESMERTONEQUATION_H_ */
