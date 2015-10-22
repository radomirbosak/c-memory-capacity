#include <stdbool.h>
#include <stdlib.h>
#include <math.h>	
#include <time.h>

#include "randomvar.h"

double unif01() {
	return (double)rand() / (double)RAND_MAX ;
}

double unif (double low, double high) {
	return low + unif01() * (high - low) ;
}

void initrand() {
	srand ( time(NULL) );
}

double saved_normal;
bool b_saved = false;
double normal01() {
	if (b_saved) {
		b_saved = false;
		return saved_normal;
	} else {
		double u = unif01();
		double v = unif01();
		double mid = sqrt(-2 * log(u));
		saved_normal = mid * cos(2 * M_PI * v);
		b_saved = true;
		return mid * sin(2 * M_PI * v);
	}
	
}

double normal(double loc, double std) {
	return loc + normal01() * std;
}