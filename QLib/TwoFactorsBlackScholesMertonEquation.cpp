/*
 * TwoFactorsBlackScholesMertonEquation.cpp
 *
 *  Created on: May 14, 2012
 *      Author: Administrator
 */

#include "TwoFactorsBlackScholesMertonEquation.h"

TwoFactorsBlackScholesMertonEquation::TwoFactorsBlackScholesMertonEquation
(double r, double vol1, double vol2, double correl)
:r(r),vol1(vol1),vol2(vol2),correl(correl)
{
	// TODO Auto-generated constructor stub
	tmp1 = 0.5*vol1*vol1;
	tmp2 = 0.5*vol2*vol2;
	tmp3 = correl*vol1*vol2;
}

TwoFactorsBlackScholesMertonEquation::~TwoFactorsBlackScholesMertonEquation()
{
	// TODO Auto-generated destructor stub
}

