#include "seidel.h"
#include "functions.h"
#include <vector>

std::pair<std::vector<double>, int>  seidel(std::vector<std::vector<double>> a, std::vector<double> b, double e)
{
	//метод Зейделя
	std::vector<std::vector<double>> c;
	int n = b.size(), count = 0;
	std::vector<double> v(n, 0);

	if (!positive(a))
	{
		std::vector<std::vector<double>> at(n, v);

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				at[i][j] = a[j][i];

		a = mult_matr(at, a);
		b = mult(at, b);
	}

	std::vector<double> vec, d;

	for (int i = 0; i < n; ++i)
	{
		c.push_back(vec);
		d.push_back(b[i] / a[i][i]);
		for (int j = 0; j < n; ++j)
		{
			if (i == j)
				c[i].push_back(0);
			else
				c[i].push_back((-a[i][j]) / a[i][i]);
		}
	}

	std::vector<double> x = d;

	do
	{
		for (int i = 0; i != n; ++i)
		{
			x[i] = 0;
			for (int j = 0; j != n; ++j)
			{
				x[i] += c[i][j] * x[j];
			}
			x[i] += d[i];
		}

		++count;
	} while (criteria_seidel(a, x, b) > e);

	return std::make_pair(x, count);
}