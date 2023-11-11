#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <vector>
#include <fstream>
void show(std::vector<double> x, std::ofstream& file_out);

void show(std::pair<std::vector<double>, int> x, std::ofstream& file_out);

double det(std::vector<std::vector<double>> A);

bool positive(std::vector<std::vector<double>> A);

double sqr(double x);

std::vector<double> mult(std::vector<std::vector<double>> a, std::vector<double> b);

std::vector<double> mult_num(double a, std::vector<double>& z);

std::vector<std::vector<double>> mult_num(double a, std::vector<std::vector<double>>& w);

std::vector<double> sum(std::vector<double> a, std::vector<double> b);

std::vector<std::vector<double>> sum(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);

double norma_vec(std::vector<double> vec);

double norma_matrix(std::vector<std::vector<double>> a);

bool check(std::vector<std::vector<double>> a);

double norma_matrix2(std::vector<std::vector<double>> a);

double norma_matrix3(std::vector<std::vector<double>> a);

std::vector<std::vector<double>> mult_matr(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);

std::vector<double> find_y(std::vector<std::vector <double>> L, std::vector <double> b);

std::vector<double> find_x(std::vector<std::vector <double>> U, std::vector <double> y);

double criteria_seidel(std::vector<std::vector<double>> a, std::vector<double> x, std::vector<double> b);

double criteria_iteration(std::vector<std::vector<double>> a, std::vector<double> b, std::vector<std::vector<double>> B, std::vector<double> x2, std::vector<double> x1, bool more_than_one);

void LU_(std::vector<std::vector<double>> A, std::vector<std::vector<double>>& U, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& P);

#endif
