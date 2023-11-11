#include "functions.h"
#include <vector>
#include <algorithm>
#include <fstream>

double det(std::vector<std::vector<double>> A)
{
	// вычисление определителя матрицы с помощью LU-разложения
	int n = A.size(), det = 0;
	std::vector<double> v(n, 0);
	std::vector<std::vector<double>> P(n, v), U(n, v), L(n, v);

	LU_(A, U, L, P);

	for (int i = 0; i < n; ++i)
		det += U[i][i];

	return det;
}

bool positive(std::vector<std::vector<double>> A)
{
	// проверка на положительную определенность матрицы
	std::vector<std::vector<double>> B;
	std::vector<double> v;
	int size = A.size();

	for (int n = 1; n <= size; ++n)
	{
		for (int i = 0; i < n; ++i)
		{
			B.push_back(v);
			for (int j = 0; j < n; ++j)
			{
				B[i].push_back(A[i][j]);
			}
		}
		if (det(B) <= 0) 
			return false;
		B.clear();
	}

	return true;
}

double sqr(double x)
{
	// вычисление корня
	double a = x, b, c;
	int k = 0;

	do
	{
		b = (a + x / a) / 2;
		c = a;
		a = b;
		++k;
	} while (abs(a - c) > 0.0000000000000001);

	return a;
}

std::vector<double> mult(std::vector<std::vector<double>> a, std::vector<double> b)
{
	// умножение матрицы на вектор 
	int n = b.size();
	std::vector<double> vec(n, 0);

	for (int i = 0; i != n; ++i)
		for (int j = 0; j != n; ++j)
			vec[i] += a[i][j] * b[j];

	return vec;
}

std::vector<double> mult_num(double a, std::vector<double>& z)
{
	// умножение вектора на число
	int n = z.size();
	for (int i = 0; i < n; ++i)
		z[i] *= a;
	return z;
}

std::vector<std::vector<double>> mult_num(double a, std::vector<std::vector<double>>& w)
{
	// умножение матрицы на число 
	int n = w.size();

	for (int i = 0; i < n; ++i)
		for (int j = 0; i < n; ++j)
			w[i][j] *= a;

	return w;
}

std::vector<double> sum(std::vector<double> a, std::vector<double> b)
{
	// сумма векторов
	int n = a.size();
	for (int i = 0; i != n; ++i)
		a[i] += b[i];
	return a;
}

std::vector<std::vector<double>> sum(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b)
{
	// сумма матриц
	int n = a.size();
	for (int i = 0; i != n; ++i)
		for (int j = 0; j != n; ++j)
			a[i][j] += b[i][j];
	return a;
}

double norma_vec(std::vector<double> vec)
{
	// норма вектора
	double norm = 0;
	int n = vec.size();

	for (int i = 0; i != n; ++i)
		norm += pow(vec[i], 2);

	if (norm == 0) 
		return 0;
	norm = sqr(norm);
	return norm;
}

double norma_matrix(std::vector<std::vector<double>> a)
{
	// норма матрицы
	double max = 0, temp = 0;
	int n = a.size();

	for (int i = 0; i != n; ++i)
	{
		for (int j = 0; j != n; ++j)
		{
			temp += abs(a[i][j]);
		}
		if (temp > max) 
			max = temp;
		temp = 0;
	}
	return max;
}

bool check(std::vector<std::vector<double>> a)
{
	int n = a.size();
	for (int i = 0; i != n; ++i)
	{
		for (int j = 0; j != n; ++j)
		{
			if (i > j && abs(a[i][j]) > 0.0000000000000001)
				return false;
		}
	}
	return true;
}

double norma_matrix2(std::vector<std::vector<double>> a)
{
	// спектральная норма матрицы
	int n = a.size();
	std::vector<double> v(n, 0);
	std::vector<std::vector<double>> P(n, v), U(n, v), L(n, v), at(n, v);

	for (int i = 0; i != n; ++i)
		for (int j = 0; j != n; ++j)
			at[i][j] = a[j][i];

	a = mult_matr(at, a);

	while (!check(a))
	{
		LU_(a, U, L, P);
		a = mult_matr(U, L);
	}

	double max = 0;
	for (int i = 0; i != n; ++i)
	{
		if (abs(a[i][i]) > max)
			max = sqrt(a[i][i]);
	}
	return max;
}

double norma_matrix3(std::vector<std::vector<double>> a)
{
	// норма матрицы
	double max = 0, temp = 0;
	int n = a.size();

	for (int j = 0; j != n; ++j)
	{
		for (int i = 0; i != n; ++i)
		{
			temp += abs(a[i][j]);
		}

		if (temp > max)
			max = temp;
		temp = 0;
	}
	return max;
}

std::vector<std::vector<double>> mult_matr(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
	// умножение матриц
	int n = A.size();
	std::vector<double> v(n, 0);
	std::vector<std::vector<double>> R(n, v);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
	return R;
}

std::vector<double> find_y(std::vector<std::vector <double>> L, std::vector <double> b)
{
	// нахождение вектора y в методе LU-разложения
	std::vector<double> y;
	int n = b.size();
	for (int i = 0; i != n; ++i)
	{
		y.push_back(b[i]);
		for (int k = 0; k < i; ++k)
		{
			y[i] -= L[i][k] * y[k];
		}
	}
	return y;
}

std::vector<double> find_x(std::vector<std::vector <double>> U, std::vector <double> y)
{
	// нахождение вектора x в методе LU - разложения и QR-разложения
	int n = y.size();
	std::vector<double> x(n, 0);

	for (int i = n - 1; i > -1; --i)
	{
		x[i] = y[i];
		for (int k = i + 1; k < n; ++k)
		{
			x[i] -= U[i][k] * x[k];

		}
		x[i] /= U[i][i];
	}
	return x;
}

double criteria_seidel(std::vector<std::vector<double>> a, std::vector<double> x, std::vector<double> b)
{
	// критерий остановки в методе Зейделя
	const int n = b.size();
	double norm = 0;
	std::vector<double> vec = mult(a, x);

	for (int i = 0; i != n; ++i)
		vec[i] -= b[i];

	return norma_vec(vec);
}

double criteria_iteration(std::vector<std::vector<double>> a, std::vector<double> b, std::vector<std::vector<double>> B, std::vector<double> x2, std::vector<double> x1, bool more_than_one)
{
	// критерий остановки в методе простой итерации
	if (more_than_one)
	{
		//если норма матрицы B больше единицы, то используется критерий ||Ax - b|| < E
		double normB = norma_matrix(B);
		std::vector<double> x;
		int n = x2.size();

		for (int i = 0; i != n; ++i)
			x.push_back(x2[i] - x1[i]);

		return (normB / (1 - normB)) * norma_vec(x);
	}

	else
	{
		//иначе критерий Зейделя
		return criteria_seidel(a, x2, b);
	}
}

void LU_(std::vector<std::vector<double>> A, std::vector<std::vector<double>>& U, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& P)
{
	// LU - разложение
	int n = A.size(), max = 0;

	for (int i = 0; i < n; ++i)
		P[i][i] = 1;

	for (int g = 0; g < n - 1; ++g)
	{
		for (int j = g; j < n; ++j)
		{
			if (abs(A[j][g]) > max)
				max = abs(A[j][g]);
		}

		for (int j = g; j < n; ++j)
		{
			if (abs(A[j][g]) == max)
			{
				max = 0;
				swap(A[j], A[g]);
				swap(P[j], P[g]);
				break;
			}
		}

		int i = g;

		for (int j = i + 1; j < n; ++j)
		{
			A[j][i] = A[j][i] / A[i][i];
			for (int k = i + 1; k < n; ++k)
			{
				A[j][k] = A[j][k] - A[j][i] * A[i][k];
			}
		}
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i <= j)
				U[i][j] = A[i][j];
			else if (i > j)
				L[i][j] = A[i][j];
			if (i == j)
				L[i][j] = 1;
		}
	}
}

void show(std::vector<double> x, std::ofstream& file_out)
{
	// вывод вектора на экран
	for (auto& c : x)
		file_out << c << " ";
}

void show(std::pair<std::vector<double>, int> x, std::ofstream& file_out)
{
	// вывод вектора на экран
	for (auto& c : x.first)
		file_out << c << " ";
}
