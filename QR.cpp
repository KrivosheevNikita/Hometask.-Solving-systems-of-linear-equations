#include "QR.h"
#include "functions.h"
#include <vector>

std::vector<double> QR(std::vector<std::vector<double>> R, std::vector<double> b)
{
	// решение методом QR
	int n = b.size();
	double a;
	std::vector<double> x(n, 0), y, v, z;
	std::vector<std::vector<double>> Q(n, x), Qn, Rn, temp;

	for (int i = 0; i < n; ++i)
		Q[i][i] = 1;

	Qn = Q;

	for (int k = 0; k < n - 1; ++k)
	{
		y.clear();

		for (int i = k; i < n; ++i)
			y.push_back(R[i][k]);

		a = norma_vec(y);
		z.clear();
		z.push_back(1);

		for (int i = 1; i < n - k; ++i)
			z.push_back(0);

		std::vector<double> w;
		mult_num(-a, z);
		w = sum(y, z);
		w = mult_num(1 / norma_vec(w), w); // wi = (yi - ai*z)/ ||(yi - ai*zi)||
		Rn.clear();
		temp.clear();

		for (int i = 0; i < n - k; ++i)
		{
			Rn.push_back(v);
			temp.push_back(v);
		}

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				if (i >= k && j >= k)
				{
					Q[i][j] = (-2 * w[i - k] * w[j - k]); // Qi = E - 2*wi*wti

					if (i == j)
						Q[i][i] += 1;

					temp[i - k].push_back(Q[i][j]);
					Rn[i - k].push_back(R[i][j]);
				}

				else
				{
					if (i == j)
						Q[i][j] = 1;
					else
						Q[i][j] = 0;
				}
			}
		}

		Qn = mult_matr(Qn, Q); // преобразование матрицы Q
		Rn = mult_matr(temp, Rn); // преобразование матрицы R

		for (int i = k; i < n; ++i)
		{
			for (int j = k; j < n; ++j)
			{
				R[i][j] = Rn[i - k][j - k];
			}
		}
	}

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			Q[i][j] = Qn[j][i];

	y = mult(Q, b);
	x = find_x(R, y);
	return x;
}
