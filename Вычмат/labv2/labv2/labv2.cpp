#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "cmath.h"
#include "spline.c"


const int COUNT_OF_POINT = 8;
const int END = 0;
const double SLOPE = 0;

double func(double x0[], double f0[], double b[], double c[], double d[], int& last, double& x);

int main()
{
	std::vector<double> x0 = {0.0, 0.2, 0.5, 0.7, 1.0, 1.3, 1.7, 2.0};
	std::vector<double> f0 = {1.0, 1.1487, 1.4142, 1.6245, 2.0, 2.4623, 3.2489, 4.0};
	std::vector<double> b(COUNT_OF_POINT);
	std::vector<double> c(COUNT_OF_POINT);
	std::vector<double> d(COUNT_OF_POINT);
	int flag = 0;

	spline(COUNT_OF_POINT, END, END, SLOPE, SLOPE, x0.data(), f0.data(), b.data(), c.data(), d.data(), &flag);

	double begin = 0;
	double end = 2;

	double eps = 1e-8;

	double dis = 0;
	while (end - begin > eps)
	{
		double center = (begin + end) / 2;
		int last = 0;

		double f_end = func(x0.data(), f0.data(), b.data(), c.data(), d.data(), last, end);
		double f_center= func(x0.data(), f0.data(), b.data(), c.data(), d.data(), last, center);

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

	double x1 = begin;
	double residual = dis;

	std::cout << "Result: " << std::scientific << x1 << " Residual: " << std::scientific << residual << '\n';
}

double func(double x0[], double f0[], double b[], double c[], double d[], int& last, double& x)
{
	return seval(COUNT_OF_POINT, x, x0, f0, b, c, d, &last) + 5 * x - 3;
}
