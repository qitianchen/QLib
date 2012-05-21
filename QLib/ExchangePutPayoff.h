/*
 * ExchangePutPayoff.h
 *
 *  Created on: May 15, 2012
 *      Author: Administrator
 */

#ifndef EXCHANGEPUTPAYOFF_H_
#define EXCHANGEPUTPAYOFF_H_

class ExchangePutPayoff
{
public:
	ExchangePutPayoff();
	virtual ~ExchangePutPayoff();
	inline double operator () (double S1, double S2) const
	{
		if ( S2 > S1 )
		{
			return S2-S1;
		}

		return 0;
	}
};

#endif /* EXCHANGEPUTPAYOFF_H_ */
