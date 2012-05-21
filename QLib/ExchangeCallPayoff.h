/*
 * ExchangeCallPayoff.h
 *
 *  Created on: May 15, 2012
 *      Author: Administrator
 */

#ifndef EXCHANGECALLPAYOFF_H_
#define EXCHANGECALLPAYOFF_H_

class ExchangeCallPayoff
{
public:
	ExchangeCallPayoff();
	virtual ~ExchangeCallPayoff();
	inline double operator () (double S1, double S2)
	{
		if ( S1 > S2 )
		{
			return S1-S2;
		}

		return 0;
	}
};

#endif /* EXCHANGECALLPAYOFF_H_ */
