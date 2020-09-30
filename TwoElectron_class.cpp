#include "project2.hpp"
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace std;
using namespace arma;


void TwoElectron::Initialize(int N, double rho_max, double omega, bool test)
{
	Solver_class::Initialize(N, test);

	m_rho_max = rho_max;
	m_omega = omega;
	m_h = m_rho_max/(m_N-1.0);

	m_Diag = 2.0/(m_h*m_h);
	m_nDiag = -1.0/(m_h*m_h);

	m_A = Dmatrix(m_Nmat, m_Nmat);
	m_R = Dmatrix(m_Nmat, m_Nmat);

	for (int i = 0; i < m_Nmat; i++) {
		for (int j = 0; j < m_Nmat; j++) {
			if (i == j) {
				m_A[i][j] = m_Diag + m_omega*m_omega*(i+1.0)*(i+1.0)*m_h*m_h + 1.0/((i+1.0)*m_h);
				m_R[i][j] = 1.0;
			}
			else {
				m_R[i][j] = 0.0;
				if (i+1==j or i-1==j)
					m_A[i][j] = m_nDiag;
				else
					m_A[i][j] = 0;
			}
		}
	}

	if (m_test)
		EigenCompare();
}

void TwoElectron::Jacobi(double tolerance, int maxiter)
{
	Solver_class::Jacobi(tolerance, maxiter);

	if (m_test)
		Jacobi_TEST();
}
