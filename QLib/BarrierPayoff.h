/*
 * BarrierPayoff.h
 *
 *  Created on: May 14, 2012
 *      Author: Administrator
 */

#ifndef BARRIERPAYOFF_H_
#define BARRIERPAYOFF_H_

class BarrierPayoff
{
public:
	BarrierPayoff(double K, double K1, double K2);
	virtual ~BarrierPayoff();
	inline double operator () (double S) const
	{
		if ( S < K )
		{
			return 0;
		}

		if ( S > K2 )
		{
			return 0;
		}

		return S - K;
	}

private:
	double K, K1, K2;
};

#endif /* BARRIERPAYOFF_H_ */
