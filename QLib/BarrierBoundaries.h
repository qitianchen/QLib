/*
 * BarrierBoundaries.h
 *
 *  Created on: May 14, 2012
 *      Author: Administrator
 */

#ifndef BARRIERBOUNDARIES_H_
#define BARRIERBOUNDARIES_H_

class BarrierBoundaries
{
public:
	BarrierBoundaries(double K1, double K2);
	virtual ~BarrierBoundaries();
	inline double lb(double t, double S) const
	{
		return 0;
	}
	inline double ub(double t, double S) const
	{
		return 0;
	}

private:
	double K1, K2;
};

#endif /* BARRIERBOUNDARIES_H_ */
