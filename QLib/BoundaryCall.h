/*
 * BoundaryCall.h
 *
 *  Created on: May 14, 2012
 *      Author: Administrator
 */

#ifndef BOUNDARYCALL_H_
#define BOUNDARYCALL_H_

#include <iostream>
#include <cmath>

class BoundaryCall
{
public:
	BoundaryCall(double K, double r);
	virtual ~BoundaryCall();
	inline double lb(double t, double minS) const
	{
		return 0;
	}
	inline double ub(double t, double maxS) const
	{
		return maxS-K*std::exp(-r*t);
	}

private:
	double K;
	double r;
};

#endif /* BOUNDARYCALL_H_ */
