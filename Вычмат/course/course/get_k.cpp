#include "get_k.hpp"
#include <cmath>
#include "cmath.h"
#include "quanc8.c"

double get_k()
{
	return 1235.802 * get_integral();
}

double func(double x)
{
	return std::sqrt((1 - 0.25 * x * x) / (1 - x * x));
}

double get_integral()
{
	double a = 0, b = 0.5;
	double abserr = 10e-3, relerr = 10e-3;
	double result = 0, errest = 0;
	int nofunr = 0;
	double posnr = 0;
	int flag = 0;

	quanc8(func, a, b, abserr, relerr, &result, &errest, &nofunr, &posnr, &flag);

	return result;
}
