#include "seidel.h"
#include "simple_iteration.h"
#include "LU.h"
#include "QR.h"
#include "functions.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>


int main() 
{
	// вывод результатов тестов в таблицы
	std::vector<double> v;
	const int N = 8;
	std::vector<std::vector<std::vector<double>>> tests_a; // матрицы A для тестов 0-4
	std::vector<std::vector<double>> tests_b, // матрицы B для тестов 0-4
		tests_x, // точные решения, полученные с помощью внешних сервисов, для тестов 0-4
		a; // матрица A для теста 5
	std::vector<double> b, // матрица B для теста 5
		x0; // точные решения, полученные с помощью внешних сервисов, для теста 5
	tests_a.push_back({{0, 2, 3},
						{1, 2, 4},
						{4, 5, 6} });
	tests_b.push_back({ 13, 17, 32 });
	tests_x.push_back({ 1, 2, 3 });

	tests_a.push_back({{N + 2, 1, 1},
						{1, N+4, 1},
						{1, 1, N+6} });
	tests_b.push_back({ N+4, N+6, N+8 });
	tests_x.push_back({ 1, 1, 1 });

	tests_a.push_back({ {-(N + 2), 1, 1},
						{1, -(N + 4), 1},
						{1, 1, -(N + 6)} });
	tests_b.push_back({ -(N + 4), -(N + 6), -(N + 8) });
	tests_x.push_back({ double(1211)/821, double(1151)/821, double(1107)/821 });

	tests_a.push_back({ {-(N + 2), N + 3, N + 4},
						{N + 5, -(N + 4), N + 1},
						{N + 4, N + 5, -(N + 6)} });
	tests_b.push_back({ (N + 4), (N + 6), (N + 8) });
	tests_x.push_back({ double(2210) / 1609, double(1840) / 1609, double(1764) / 1609 });

	tests_a.push_back({ {N + 2, N + 1, N + 1},
						{N + 1, N + 4, N + 1},
						{N + 1, N + 1, N + 6} });
	tests_b.push_back({ (N + 4), (N + 6), (N + 8) });
	tests_x.push_back({ double(-3) / 37, double(71) / 111, double(29) / 37 });

	std::ofstream file_out1, file_out2;
	file_out1.open("Excel_table3.csv"); // результаты тестов 0-4
	file_out2.open("Excel_table4.csv"); // результаты теста 5

	// вывод результатов тестов 0-4 в таблицу
	file_out1 << "№ теста;x0;e;Метод Зейделя;;;;Метод простой итерации;;;;Метод Гаусса;;;Метод Хаусхолдера" << std::endl;
	file_out1 << ";;;x;d;k;;x;d;k;;x;d;;x;d" << std::endl;
	for (int i = 0; i < 5; ++i) 
	{
		file_out1 << i << ";";
		show(tests_x[i], file_out1);

		for (double e = 0.01; e > 0.000001; e*=0.1) 
		{
			file_out1 <<";" << e << ";";
			auto x = seidel(tests_a[i], tests_b[i], e);
			show(x.first, file_out1);

			file_out1 << ";" << norma_vec(sum(mult_num(-1, x.first), tests_x[i])) << ";" << x.second;
			file_out1 << ";;";
			x = simple_iteration(tests_a[i], tests_b[i], e);
			show(x.first, file_out1);

			file_out1 << ";" << norma_vec(sum(mult_num(-1, x.first), tests_x[i])) << ";" << x.second;

			if (e == 0.01) 
			{
				file_out1 << ";;";
				auto x = LU(tests_a[i], tests_b[i]);
				show(x, file_out1);

				file_out1 << ";" << norma_vec(sum(mult_num(-1, x), tests_x[i])) << ";;";
				x = QR(tests_a[i], tests_b[i]);
				show(x, file_out1);

				file_out1 << ";" << norma_vec(sum(mult_num(-1, x), tests_x[i]));
			}

			file_out1 << std::endl;
			file_out1 << ";";
		}
		file_out1 << std::endl << std::endl;
	}

	file_out1.close();
	system("Excel_table3.csv");
	// вывод результатов теста 5 в таблицу
	file_out2 << "5 тест" << std::endl;
	file_out2 << "n;E;x0;e;Метод Зейделя;;;;Метод простой итерации;;;;Метод Гаусса;;;Метод Хаусхолдера" << std::endl;
	file_out2 << ";;;;x;d;k;;x;d;k;;x;d;;x;d" << std::endl;

	for (int n = 4; n < 11; ++n) 
	{
		file_out2 << n;

		for (double E = 0.001; E > 0.0000001; E *= 0.001) 
		{
			file_out2 << ";" << E << ";";
			a.clear();
			b.clear();
			x0.clear();

			for (int i = 0; i < n; ++i) 
			{
				a.push_back(v);
				if (i != n - 1) 
					b.push_back(-1);
				else 
					b.push_back(1);
				if (i != n - 1) 
					x0.push_back(0);
				else if (E == 0.001) 
					x0.push_back(double(1000) / 1001);
				else  x0.push_back(double(1000000) / 1000001);

				for (int j = 0; j < n; ++j) 
				{
					if (i == j) 
						a[i].push_back(1 + E);
					else if (i > j) 
						a[i].push_back(E);
					else if (i < j) 
						a[i].push_back(-1 - E);
				}
			}

			for (double e = 0.01; e > 0.000001; e *= 0.1) 
			{
				if (e == 0.01) 
				{
					show(x0, file_out2);
					file_out2 << ";";
				}

				file_out2 << e << ";";
				auto x = seidel(a, b, e);
				show(x, file_out2);

				file_out2 << ";" << norma_vec(sum(x.first, mult_num(-1, x0)));
				mult_num(-1, x0);
				file_out2 << ";" << x.second << ";;";
				x = simple_iteration(a, b, e);
				show(x, file_out2);

				file_out2 << ";" << norma_vec(sum(x.first, mult_num(-1, x0)));
				mult_num(-1, x0);
				file_out2 << ";" << x.second << ";;";
				if (e == 0.01)
				{
					auto x = LU(a, b);
					show(x, file_out2);

					file_out2 << ";" << norma_vec(sum(x, mult_num(-1, x0))) << ";;";
					mult_num(-1, x0);
					x = QR(a, b);
					show(x, file_out2);

					file_out2 << ";" << norma_vec(sum(x, mult_num(-1, x0))) << ";;";
					mult_num(-1, x0);
				}
				file_out2 << std::endl;
				if (e > 0.000001) 
					file_out2 << ";;;";
			}
			file_out2 << std::endl;
		}
	}

	file_out2.close();
	system("Excel_table4.csv");
    return 0;
	}
