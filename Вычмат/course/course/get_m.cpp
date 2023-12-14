#include "get_m.hpp"
#include <cmath>

double get_m()
{
	double x = bisection(func, 0, 2, 10e-4);
	return x * 25.04499;
}

double func(double& x)
{
	return std::pow(2, x) - 2 * x * x - 1;
}

double bisection(double (*func)(double& x), double begin, double end, double eps)
{
	double dis = 0;
	while (end - begin > eps)
	{
		double center = (begin + end) / 2;
		int last = 0;

		double f_end = func(end);
		double f_center = func(center);
		if (f_end * f_center < 0)
		{
			begin = center;
		}
		else
		{
			end = center;
		}
		dis = f_center;
	}

	return begin;
}