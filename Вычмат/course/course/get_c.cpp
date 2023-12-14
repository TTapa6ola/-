#include "get_c.hpp"
#include "cmath.h"
#include "decomp.c"

Matrix a = {46, 42, 24, 42, 49, 18, 24, 18, 16};
Matrix b = {282, 229, 178};

const int size = 3;
double cond = 0;

Matrix get_c()
{
	Matrix x(size);
	Matrix inv_a = create_inverse_matrix(a);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			x[i] += inv_a[i * size + j] * b[j];
		}
	}

	return x;
}

Matrix create_inverse_matrix(Matrix m)
{
	Matrix matrix(size * size);

	int pivot[size];
	int flag = 0;

	decomp(size, size, m.data(), &cond, pivot, &flag);

	Matrix e = create_identity_matrix();

	for (int i = 0; i < size; i++)
	{
		Matrix b(size);

		for (int j = 0; j < size; j++)
		{
			b[j] = e[j * size + i];
		}

		solve(size, size, m.data(), b.data(), pivot);

		for (int k = 0; k < size; k++)
		{
			matrix[k * size + i] = b[k];
		}
	}

	return matrix;
}

Matrix create_identity_matrix()
{
	Matrix m(size * size);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			m[i * size + j] = (i == j) ? 1 : 0;
		}
	}

	return m;
}
