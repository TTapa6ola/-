#include <iostream>
#include <vector>
#include "get_c.hpp"
#include "get_m.hpp"
#include "get_k.hpp"
#include "rkf45.c"

Matrix R = {80, 160, 240};
Matrix C;
double M, K;
struct pair
{
	double r0, c;
};

std::vector<pair> get_pairs_of_r_c(Matrix C, Matrix R);
void do_rkf(int (*f)(int n, double t, double x[], double xp[]));
double x0(double t);
double dx0dt(double t);
double r1(double t, double x);
int func1(int n, double t, double x[], double dxdt[]);
double r2(double t, double x);
int func2(int n, double t, double x[], double dxdt[]);
double r3(double t, double x);
int func3(int n, double t, double x[], double dxdt[]);


int main()
{
	M = get_m();
	K = get_k();
	C = get_c();
	std::vector<pair> pairs_r_c = get_pairs_of_r_c(get_c(), R);

	std::cout << "M = " << M << ", " << "K = " << K << '\n';
	for (int i = 0; i < 3; i++)
	{
		std::cout << "r" << i + 1 << " = " << pairs_r_c.at(i).r0 << ", c" << i + 1 << " = " << pairs_r_c.at(i).c << '\n';
	}


	std::cout << "r1, c1\n";
	do_rkf(func1);
	std::cout << "r2, c2\n";
	do_rkf(func2);
	std::cout << "r3, c3\n";
	do_rkf(func3);
}

void do_rkf(int (*f)(int n, double t, double x[], double xp[]))
{
	Matrix x_start = {0, 0};
	Matrix dxdt(2);
	double t_init = 0, out = 0, re = 1e-5, ae = 1e-5, h = 0;
	int flag = 1, fail = 1, nfe = 0, maxnfe = 10e7;
	rkfinit(2, &fail);

	for (double t = 0; t <= 4.001; t += 0.04)
	{
		rkf45(f, 2, x_start.data(), dxdt.data(), &t_init, t, &re, ae, &h, &nfe, maxnfe, &flag);
		std::cout << t << '\t' << x_start[0] << ' ' << x_start[1] << '\n';
	}

	//rkfend;
}

std::vector<pair> get_pairs_of_r_c(Matrix C, Matrix R)
{
	std::vector<pair> pairs;

	for (int i = 0; i < 3; i++)
	{
		pairs.push_back({R.at(i), C.at(i)});
	}

	return pairs;
}

double x0(double t)
{
	return 2 * (1 - cos(7 * t));
}

double dx0dt(double t)
{
	return 14 * sin(t);
}

double r1(double t, double x)
{
	return R[0] * (1 + C[0] * abs(x - dx0dt(t)));
}

int func1(int n, double t, double x[], double dxdt[])
{
	dxdt[0] = (K * (x[1] - x0(t)) + r1(t, x[0]) * (x[0] - dx0dt(t))) / (-M);
	dxdt[1] = x[0];
	return 0;
}

double r2(double t, double x)
{
	return R[1] * (1 + C[1] * abs(x - dx0dt(t)));
}

int func2(int n, double t, double x[], double dxdt[])
{
	dxdt[0] = (K * (x[1] - x0(t)) + r2(t, x[0]) * (x[0] - dx0dt(t))) / (-M);
	dxdt[1] = x[0];
	return 0;
}

double r3(double t, double x)
{
	return R[2] * (1 + C[2] * abs(x - dx0dt(t)));
}

int func3(int n, double t, double x[], double dxdt[])
{
	dxdt[0] = (K * (x[1] - x0(t)) + r3(t, x[0]) * (x[0] - dx0dt(t))) / (-M);
	dxdt[1] = x[0];
	return 0;
}
