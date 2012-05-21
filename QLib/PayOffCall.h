/*
 * PayOffCall.h
 *
 *  Created on: May 14, 2012
 *      Author: Administrator
 */

#ifndef PAYOFFCALL_H_
#define PAYOFFCALL_H_

class PayOffCall
{
public:
	PayOffCall(double K);
	virtual ~PayOffCall();
	inline double operator () (double S) const
	{
		if ( S < K )
		{
			return 0;
		}

		return S-K;
	}

private:
	double K;
};

#endif /* PAYOFFCALL_H_ */
