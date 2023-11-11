#include "simple_iteration.h"
#include "functions.h"
#include <vector>

std::pair<std::vector<double>, int> simple_iteration(std::vector<std::vector<double>> a, std::vector<double> b, double e)
{
	// метод простой итерации
	int n = b.size(), count = 0;
	std::vector<double> v(n, 0);
	bool more_than_one = false;
	std::vector<std::vector<double>> B;
	std::vector<double> vec, c;
	double t = (1 / (norma_matrix(a)));

	for (int i = 0; i != n; ++i)
	{
		B.push_back(vec);
		for (int j = 0; j != n; ++j)
		{
			if (i == j)
				B[i].push_back(1 - (t * a[i][j]));
			else
				B[i].push_back(-a[i][j] * t);
		}
	}

	if (norma_matrix(B) > 1 && norma_matrix2(B) > 1 && norma_matrix3(B) > 1)
	{
		//если норма больше единицы, то домножается на транспонированную матрицу
		std::vector<std::vector<double>> at(n, v);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				at[i][j] = a[j][i];

		a = mult_matr(at, a);
		b = mult(at, b);
		B.clear();
		t = (1 / (norma_matrix(a)));

		for (int i = 0; i != n; ++i)
		{
			B.push_back(vec);
			for (int j = 0; j != n; ++j)
			{
				if (i == j)
					B[i].push_back(1 - (t * a[i][j]));
				else
					B[i].push_back(-a[i][j] * t);
			}
		}
	}

	if (norma_matrix(B) > 1 && norma_matrix2(B) > 1 && norma_matrix3(B) > 1)
		more_than_one = true; // если норма больше единицы, то используется критерий остановки из метода Зейделя

	for (int i = 0; i != n; ++i)
		c.push_back(t * b[i]);

	std::vector<double> x = c, temp;

	do
	{
		temp = x;
		x = sum(mult(B, x), c);
		++count;
	} while (criteria_iteration(a, b, B, x, temp, more_than_one) > e);

	return std::make_pair(x, count);
}
