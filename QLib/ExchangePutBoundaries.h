/*
 * ExchangePutBoundaries.h
 *
 *  Created on: May 15, 2012
 *      Author: Administrator
 */

#ifndef EXCHANGEPUTBOUNDARIES_H_
#define EXCHANGEPUTBOUNDARIES_H_

class ExchangePutBoundaries
{
public:
	ExchangePutBoundaries(double S1Max, double S2Max);
	virtual ~ExchangePutBoundaries();
	inline double lbx(double t, double minS1, double S2) const
	{
		return S2;
	}
	inline double ubx(double t, double maxS1, double S2) const
	{
		return 0;
	}
	inline double lby(double t, double S1, double minS2) const
	{
		return 0;
	}
	inline double uby(double t, double S1, double maxS2) const
	{
		return S2Max-S1;
	}

private:
	double S1Max;
	double S2Max;
};

#endif /* EXCHANGEPUTBOUNDARIES_H_ */
