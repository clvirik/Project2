#pragma once
#include <string>

class Solver_class {
public:
	int m_N, m_Nmat, m_Npoint, m_iterations;
	bool m_test;
	double m_h, m_Diag, m_nDiag;
	double** m_A, ** m_R;
	double* m_eigvals, ** m_eigvecs;
	int m_k, m_l;

	double m_pi = 3.14159265358979323846;

	double* EigenVals_arma();
	double* EigenVals_anal();

	double max_OffDiag();
	void Jacobi_rotate();
	void Jacobi_sort();

	double** Dmatrix(int row, int col);
	void delete_Dmatrix(double** M, int row, int col);

	void Initialize(int N, bool test);
	void Jacobi(double tolerance, int maxiter);
};

class Basic : public Solver_class {
private:
	void EigenCompare(double eps = 1.0E-10);
	void Jacobi_TEST(double eps = 1.0E-10);

public:

	void Initialize(int N, bool test = false);
	void Jacobi(double tolerance = 1.0E-10, int maxiter = (int)1.0E+6);
};

class OneElectron : public Solver_class {
private:
	double m_rho_max;

	void EigenCompare(double eps = 1.0E-10);
	void Jacobi_TEST(double eps = 1.0E-10);

public:

	void Initialize(int N, double rho_max, bool test = false);
	void Jacobi(double tolerance = 1.0E-10, int maxiter = (int)1.0E+6);
};

class TwoElectron : public Solver_class {
private:
	double m_rho_max, m_omega;

	void EigenCompare(double eps = 1.0E-10);
	void Jacobi_TEST(double eps = 1.0E-10);

public:

	void Initialize(int N, double rho_max, double omega, bool test = false);
	void Jacobi(double tolerance = 1.0E-10, int maxiter = (int)1.0E+6);
};