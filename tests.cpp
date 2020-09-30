#include "project2.hpp"
#include <iostream>
#include <iomanip>

using namespace std;

void Basic::EigenCompare(double eps)
{
	double* eigen_anal, * eigen_arma;
	eigen_anal = new double[m_Nmat];
	eigen_arma = new double[m_Nmat];

	eigen_anal = EigenVals_anal();
	eigen_arma = EigenVals_arma();

	cout << "Comparing analytic and armadillo eigenvalues:" << endl;
	cout << "Index      Armadillo  Analytic" << endl;

	bool compare;
	int index = 0;
	for (int i = 0; i<m_Nmat; i++) {
		compare = fabs(eigen_arma[i] - eigen_anal[i]) >= eps;

		if (compare) {
			cout << "index = " << i << ": ";
			cout << eigen_arma[i] << " != " << eigen_anal[i] << endl;
			index++;
		}
	}
	delete[] eigen_anal;
	delete[] eigen_arma;

	cout << endl;
	if (not index)
		cout << "SUCCESS; THEY ARE EQUAL! epsilon = " << eps << endl << endl << endl;
	else
		cout << endl << endl;
}

void Basic::Jacobi_TEST(double eps)
{
	Jacobi_sort();

	double** vectors = Dmatrix(m_Nmat, m_Nmat);

	for (int i = 0; i<m_Nmat; i++) {
		for (int j = 0; j<m_Nmat; j++) {
			vectors[i][j] = sin((i+1.0)*(j+1.0)*m_pi/(double)m_N);
		}
	}

	cout << "Comparing Jacobi and analytic eigenvectors:" << endl;

	bool compare;
	int index = 0;
	double dot_value, analVector_len, JacobiVector_len;

	for (int i = 0; i<m_Nmat; i++) {
		dot_value = 0.0; analVector_len = 0.0; JacobiVector_len = 0.0;
		for (int j = 0; j<m_Nmat; j++) {
			dot_value += vectors[i][j]*m_eigvecs[i][j];
			analVector_len += vectors[i][j]*vectors[i][j];
			JacobiVector_len += m_eigvecs[i][j]*m_eigvecs[i][j];
		}
		compare = fabs(dot_value*dot_value - analVector_len*JacobiVector_len) >= eps;

		if (compare) {
			cout << "index = " << i << ": " << endl
				<< "Analytic     Jacobi" << endl;
			for (int j=0;j<m_Nmat;j++)
				cout << setw(10) << vectors[i][j]*m_eigvecs[i][0]/vectors[i][0] << "   " << m_eigvecs[i][j] << endl;
			index++;
		}
	}

	cout << endl;
	if (not index)
		cout << "SUCCESS; THEY ARE EQUAL! epsilon = " << eps << endl << endl << endl;
	else
		cout << endl << endl;


	for (int i = 0; i<m_Nmat; i++) {
		for (int j = 0; j<m_Nmat; j++) {
			cout << setw(10) << m_R[i][j] << "  ";
		}
		cout << endl;
	}
	cout << endl;

	delete[] m_eigvals;
	delete[] m_eigvecs;
	delete_Dmatrix(vectors, m_Nmat, m_Nmat);
}

void OneElectron::EigenCompare(double eps)
{

}

void TwoElectron::EigenCompare(double eps)
{

}

void OneElectron::Jacobi_TEST(double eps)
{

}

void TwoElectron::Jacobi_TEST(double eps)
{

}