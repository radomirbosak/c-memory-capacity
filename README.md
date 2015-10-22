Memory capacity
===============

A `C` implementation of the memory capacity function. Calculates the memory capacity of an echo-state network (a type of artificial neural network) as defined by Jeager (2001) in paper *Short term memory in echo state networks*. In this case, an instance of an echo-state network is given by the recurrent matrix W and input matrix W<sup>I</sup>.

Unfortunately, it is slower than the python implementation (which uses `numpy`), probably because of using custom vector and matrix operations.

Usage
-----
Call function `MC`

	double MC(double* W, double* WI, double sigma, double tau, int memory_max);

defined in the file `foo.c`.


Dependencies
------------
Gsl (GNU scientific library)
