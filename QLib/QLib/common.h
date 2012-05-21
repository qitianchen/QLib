/*!
 * \file common.h
 *
 * \brief Some common macros to be used across the library
 *
 *  Created on: May 1, 2012
 *      Author: Qitian Chen
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <cmath>

//! Used for comparing two double numbers
#define PRECISION 1e-10
//! Used for numeric solvers, to prevent infinite loops
#define MAX_ITERATIONS 32766

inline int round_double_to_int(double x)
{
	double ux = std::ceil(x), lx = std::floor(x);
	return (int)(((ux-x)<(x-lx))?ux:lx);
}

#endif /* COMMON_H_ */
