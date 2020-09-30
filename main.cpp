#include "project2.hpp"
#include <iostream>
#include <iomanip>
#include "time.h"
#include <fstream>

// I apologize for the quality of this project,
// but it was rather rushed as you will be able to tell

using namespace std;

int main() {
	clock_t start, finish;

	//Runs Unit Tests on the basic case
	if (false) {
		Basic my_solver;
		my_solver.Initialize(5, true);
		my_solver.Jacobi();
	}

	//Runs Unit Tests on the one electron case
	if (false) {

	}

	//Runs Unit Tests on the basic case
	if (false) {

	}

	if (false) {
		int length = 5;
		int iterations = 10;
		int Jacobi_iterations;

		double time_Jacobi, time_Arma;

		ofstream outFile;
		outFile.open("Time.txt");

		for (int N = 1; N<=length; N += 1) {
			start = clock();
			for (int i = 0; i<iterations; i++) {
				Basic my_solver;
				my_solver.Initialize(N*100);
				my_solver.EigenVals_arma();
			}
			finish = clock();
			time_Arma = (double)(finish - start)/(double)iterations;

			start = clock();
			for (int i = 0; i<iterations; i++) {
				Basic my_solver;
				my_solver.Initialize(N*100);
				my_solver.Jacobi();
				Jacobi_iterations = my_solver.m_iterations;
			}
			finish = clock();
			time_Jacobi = (double)(finish - start)/(double)iterations;

			outFile << setprecision(8) << N << " " << Jacobi_iterations << " " << time_Jacobi << " " << time_Arma << endl;
		}

		outFile.close();
	}

	if (false) {

		ofstream outFile;
		outFile.open("OneElectron.txt");
		for (int N=100; N<=300; N+=50){
			for (double rho = 2.0; rho<=6.0; rho+=0.2) {
				OneElectron my_solver;
				my_solver.Initialize(N, rho);
				my_solver.Jacobi();
				my_solver.Jacobi_sort();

				outFile << setprecision(8) << N << " " << rho << " " << fabs(my_solver.m_eigvals[0] - 3.0) << endl;
			}
		}

		outFile.close();
	}


	return 0;
}
