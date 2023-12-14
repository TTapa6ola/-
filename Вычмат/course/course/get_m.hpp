#pragma once

double get_m();
double func(double& x);
double bisection(double (*func)(double& x), double a, double b, double eps);