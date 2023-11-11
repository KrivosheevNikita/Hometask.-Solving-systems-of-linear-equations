#include "LU.h"
#include "functions.h"
#include <vector>

std::vector<double> LU(std::vector<std::vector<double>> R, std::vector<double> b)
{
	// решение методом LU
	int n = b.size();
	std::vector<double> v(n, 0);
	std::vector<std::vector<double>> P(n, v), U(n, v), L(n, v);
	LU_(R, U, L, P);
	b = mult(P, b);
	return find_x(U, find_y(L, b));
}
