QLib
====

Numerical methods for financial computing

================================================================

This project is inspired by FE8822 Numerical Methods (Finite Difference Methods), a course in Master of Science in Financial Engineering, Nanyang Technological University. The syllabus includes,
1) Non-linear equation
2) System of linear equations
3) Ordinary differential equations
4) Partial differential equations

Almost all the numerical routines are rewritten in C++ with object oriented design in mind while maintaining readability. Highlights include explicit, implicit and Crank-Nicolson methods for Black-Scholes equations. A more efficient algorithm, Operator Splitting, is implemented to address "the curse of dimensionality" in multi-asset option pricing with cross derivative term. The library is released with pricing examples of European call option, barrier option and exchange put option. By design, the library will also be able to handle practical valuation problems with yield curve and volatility surface.

The library was built on Eigen 3, 
http://eigen.tuxfamily.org/index.php?title=Main_Page
a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.


