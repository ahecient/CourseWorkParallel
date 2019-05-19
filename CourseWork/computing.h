#pragma once
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
class computing{
public:
	computing() : length(1.0f), x_points_count(30), t_points_count(2000), a(2), b(2), A(-3), B(-1),

		step_x(length / (x_points_count)), step_t(length / (t_points_count)),
		error((step_t) / (step_x * step_x)), r(a * error) {}
	~computing();

	std::vector<std::vector<long double>> approximate_parallel;
	std::vector<std::vector<long double>> precision_parallel;

	std::vector<std::vector<long double>> approximate_serial;
	std::vector<std::vector<long double>> precision_serial;

	void preparation();
	void work(const bool, std::vector<std::vector<long double>>&, std::vector<std::vector<long double>>&);
	void computing_epsilon();
	

private:
	const long double length;
	const int x_points_count;
	const int t_points_count;
	const long double step_x;
	const long double step_t;
	const double long error;
	const double long r;
	const long double A;
	const long double B;
	const long double a;
	const long double b;

	std::vector<std::vector<long double>> absolute_parallel;
	std::vector<std::vector<long double>> relative_parallel;

	void fill_vec(std::vector<std::vector<long double>>&);
	long double precision_computing(const long double x, const long double t, const long double b, const long double A, const long double B);
	long double approximate_computing(const long double w, const long double wleft, const long double wright, const long double r, const long a, const long b);
	long double maximum(std::vector<std::vector<long double>> &);
	void epsilon(std::vector<std::vector<long double>>&, std::vector<std::vector<long double>>&);
	void write();
};


