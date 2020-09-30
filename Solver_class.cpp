#include "project2.hpp"
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace std;
using namespace arma;


void Solver_class::Initialize(int N, bool test)
{
	m_N = N;
	m_Npoint = m_N+1;
	m_Nmat = m_N-1;

	m_test = test;
}

double* Solver_class::EigenVals_arma()
{
	// Getting eigenvalues using armadillo
	mat A(m_Nmat, m_Nmat);
	for (int i = 0; i < m_Nmat; i++) {
		for (int j = 0; j < m_Nmat; j++)
			A(i, j) = m_A[i][j];
	}

	vec eigval;
	mat eigvec;

	eig_sym(eigval, eigvec, A);

	double* eig = new double[m_Nmat];
	for (int i=0; i < m_Nmat; i++)
		eig[i] = eigval(i);

	return eig;
	delete[] eig;
}

double* Solver_class::EigenVals_anal()
{
	// Getting eigenvalues from analytic expression

	double* eig = new double[m_Nmat];	
	for (int i=0; i < m_Nmat; i++)
		eig[i] = m_Diag + 2*m_nDiag*cos((i+1.0)*m_pi/(double)m_N);

	return eig;
	delete[] eig;
}

void Solver_class::Jacobi(double tolerance, int maxiter)
{
	m_iterations = 0;
	while (max_OffDiag() > tolerance && m_iterations <= maxiter) {
		Jacobi_rotate(); m_iterations++;
	}
}

void Solver_class::Jacobi_sort()
{
	m_eigvals = new double[m_Nmat];
	m_eigvecs = Dmatrix(m_Nmat, m_Nmat);

	for (int i = 0; i<m_Nmat; i++)
		m_eigvals[i] = m_A[i][i];

	sort(m_eigvals, m_eigvals+m_Nmat);

	int k;
	for (int i = 0; i<m_Nmat; i++) {
		for (int j = 0; j<m_Nmat; j++) {
			if (m_eigvals[i] == m_A[j][j])
				k = j;
		}
		for (int j = 0; j<m_Nmat; j++)
			m_eigvecs[i][m_Nmat-1-j] = m_R[j][k];
	}

	delete_Dmatrix(m_A, m_Nmat, m_Nmat);
	delete_Dmatrix(m_R, m_Nmat, m_Nmat);
}

void Solver_class::Jacobi_rotate()
{
	double s, c;
	if (m_A[m_k][m_l] != 0.0) {
		double t, tau;
		tau = (m_A[m_l][m_l] - m_A[m_k][m_k])/(2*m_A[m_k][m_l]);

		if (tau >= 0) {
			t = 1.0/(tau + sqrt(1.0 + tau*tau));
		}
		else {
			t = -1.0/(-tau +sqrt(1.0 + tau*tau));
		}

		c = 1/sqrt(1+t*t);
		s = c*t;
	}
	else {
		c = 1.0;
		s = 0.0;
	}
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = m_A[m_k][m_k];
	a_ll = m_A[m_l][m_l];
	m_A[m_k][m_k] = c*c*a_kk - 2.0*c*s*m_A[m_k][m_l] + s*s*a_ll;
	m_A[m_l][m_l] = s*s*a_kk + 2.0*c*s*m_A[m_k][m_l] + c*c*a_ll;
	m_A[m_k][m_l] = 0.0;  // hard-coding non-diagonal elements by hand
	m_A[m_l][m_k] = 0.0;  // same here
	for (int i = 0; i < m_Nmat; i++) {
		if (i != m_k && i != m_l) {
			a_ik = m_A[i][m_k];
			a_il = m_A[i][m_l];
			m_A[i][m_k] = c*a_ik - s*a_il;
			m_A[m_k][i] = m_A[i][m_k];
			m_A[i][m_l] = c*a_il + s*a_ik;
			m_A[m_l][i] = m_A[i][m_l];
		}
		//  And finally the new eigenvectors
		r_ik = m_R[i][m_k];
		r_il = m_R[i][m_l];

		m_R[i][m_k] = c*r_ik - s*r_il;
		m_R[i][m_l] = c*r_il + s*r_ik;
	}
}

double Solver_class::max_OffDiag()
{
	double max = 0.0;
	for (int i = 0; i<m_Nmat; i++) {
		for (int j = i+1; j<m_Nmat; j++) {
			double aij = fabs(m_A[i][j]);
			if (aij > max) {
				max = aij; m_k = i; m_l = j;
			}
		}
	}

	return max;
}

double** Solver_class::Dmatrix(int row, int col)
{
	double** M = new double* [row];
	for (int i = 0; i < row; i++)
		M[i] = new double[col];

	return M;
}

void Solver_class::delete_Dmatrix(double** M, int row, int col)
{
	for (int i = 0; i < row; i++)
		delete[] M[i];

	delete[] M;
}