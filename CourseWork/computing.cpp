#include "computing.h"


computing::~computing(){
}

void computing::fill_vec(std::vector<std::vector<long double>> &vec){
	vec.resize(x_points_count + 1);

	for (int i = 0; i < vec.size(); ++i){
		vec.at(i).resize(t_points_count);
		std::fill(vec.at(i).begin(), vec.at(i).end(), 0);
	}
}

long double computing::precision_computing(const long double x, const long double t, const long double b, const long double A, const long double B){
	return	(exp(b * t) * pow((A * x + B), -1));
}

long double computing::approximate_computing(const long double w, const long double wleft, const long double wright, const long double r, const long double step_t, const long a, const long b){
	return (w + a * r * pow(w, -2)*((wleft - 2 * w + wright) - pow((wright - wleft), 2)/(2*w)) + step_t * b * w);
}

void computing::preparation(){
	if (r > 0.5){
		std::cerr << "ERROR. r =  " << r << std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}

	fill_vec(approximate_parallel);
	fill_vec(precision_parallel);
	fill_vec(approximate_serial);
	fill_vec(precision_serial);

	for (int i = 0; i < x_points_count + 1; ++i){
		approximate_parallel[i][0] = precision_computing(i * step_x, 0, b, A, B);
		approximate_serial[i][0] = precision_computing(i * step_x, 0, b, A, B);

	}

	for (int j = 1; j < t_points_count; ++j){
		approximate_parallel[0][j] = precision_computing(0, j * step_t, b, A, B);
		approximate_parallel[x_points_count][j] = precision_computing(x_points_count * step_x, j * step_t, b, A, B);

		approximate_serial[0][j] = precision_computing(0, j * step_t, b, A, B);
		approximate_serial[x_points_count][j] = precision_computing(x_points_count * step_x, j * step_t, b, A, B);
	}
}

void computing::work(const bool flag, std::vector<std::vector<long double>> &vec_approximate, std::vector<std::vector<long double>> &vec_precision){
	int threads_num = 0;
#pragma omp parallel if (flag)
	{
	for (int j = 1; j < t_points_count; ++j){
#pragma omp for
		for (int i = 1; i < x_points_count; ++i){
			vec_approximate[i][j] = approximate_computing(vec_approximate[i][j - 1], vec_approximate[i - 1][j - 1], vec_approximate[i + 1][j - 1], r, step_t, a, b);
			vec_precision[i][j] = precision_computing(i * step_x, j * step_t, b, A, B);
		}
	}
	threads_num = omp_get_num_threads();
	}
	std::cout << "Number of threads: " << threads_num << std::endl;
}

long double computing::maximum(std::vector<std::vector<long double>> &vec){
	long double  maximum = std::numeric_limits<long double>::min();
	for (int i = 0; i < vec.size(); ++i){
		long double temp_maximum = *std::max_element(vec.at(i).begin(), vec.at(i).end());
		if (temp_maximum > maximum){
			maximum = temp_maximum;
		}
	}
	return  (maximum);
}

void computing::epsilon(std::vector<std::vector<long double>>&vec_absolute, std::vector<std::vector<long double>>&vec_relative){
	for (int j = 1; j < t_points_count; ++j){
		for (int i = 1; i < x_points_count; ++i){
			vec_absolute[i][j] = fabs(precision_parallel[i][j] - approximate_parallel[i][j]);
			vec_relative[i][j] = vec_absolute[i][j] / precision_parallel[i][j];
		}
	}
	long double maximum_absolute = maximum(vec_absolute);
	std::cout << "Absolute: " << maximum_absolute << std::endl;
	long double maximum_relative = maximum(vec_relative);
	std::cout << "Relative: " << maximum_relative * 100 << "%" << std::endl;
}

void computing::write(){
	std::ofstream fout("results.txt");

	for (int i = 0; i < x_points_count + 1; ++i){
		for (int j = 0; j < t_points_count + 1; ++j){
			fout << i << "\t" << j << "\t" << std::setprecision(17) << std::fixed << approximate_parallel[i][j] << std::endl;
		}
	}
	std::cout << "finish writing in the file." << std::endl;
}

void computing::computing_epsilon(){
	
	fill_vec(absolute_parallel);
	fill_vec(relative_parallel);
	epsilon(absolute_parallel, relative_parallel);
	
	write();
}
