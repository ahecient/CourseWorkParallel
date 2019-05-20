#include "computing.h"
#include <chrono>
#include <conio.h>


int main(int argc, char **argv){
	computing computing;

	computing.preparation();

	auto start_ser = std::chrono::system_clock::now();
	computing.work(false, computing.approximate_serial, computing.precision_serial);
	auto end_ser = std::chrono::system_clock::now();

	auto duration_ser = std::chrono::duration_cast<std::chrono::milliseconds>(end_ser - start_ser).count();
	std::cout << "Serial : " << duration_ser << " ms" << std::endl;

	auto start_par = std::chrono::system_clock::now();
	computing.work(true, computing.approximate_parallel, computing.precision_parallel);
	auto end_par = std::chrono::system_clock::now();

	auto duration_par = std::chrono::duration_cast<std::chrono::milliseconds>(end_par - start_par).count();
	std::cout << "Parallel : " << duration_par << " ms" << std::endl;

	computing.computing_epsilon();

	_getch();
	return (0);
}