#include <math.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <complex.h>
#include <cstring>
#include <vector>

#define s_complex float _Complex
#define s_Complex(x,y) (float)x + (float)y * _Complex_I

using namespace std;

extern "C" void dsaupd_(int* ido, char* bmat, int* n, char* which,
	int* nev, double* tol, double* resid, int* ncv,
	double* v, int* ldv, int* iparam, int* ipntr,
	double* workd, double* workl, int* lworkl,
	int* info);
extern "C" void dseupd_(int* rvec, char* All, int* select, double* d,
	double* z, int* ldz, double* sigma,
	char* bmat, int* n, char* which, int* nev,
	double* tol, double* resid, int* ncv, double* v,
	int* ldv, int* iparam, int* ipntr, double* workd,
	double* workl, int* lworkl, int* ierr);
extern "C" void cnaupd_(int* ido, char* bmat, int* n, char* which,
	int* nev, double* tol, s_complex * resid, int* ncv,
	s_complex * v, int* ldv, int* iparam, int* ipntr,
	s_complex * workd, s_complex * workl, int* lworkl, double* rwork,
	int* info);
extern "C" void cneupd_(int* rvec, char* All, int* select, s_complex * d,
	s_complex * z, int* ldz, double* sigmar, double* sigmai,
	char* bmat, int* n, char* which, int* nev,
	double* tol, s_complex * resid, int* ncv, s_complex * v,
	int* ldv, int* iparam, int* ipntr, s_complex * workd,
	s_complex * workl, int* lworkl, double* rwork, int* ierr);

bool fact = false, real_size = false; //print matrix
int numbermatr = 3; //выводимая часть

double** double_A;
s_complex** complex_A;
int tlen;
void av(int n, double* in, double* out)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		out[i] = 0;
	}
	for (i = 0; i < tlen; i++)
	{
		out[(int)double_A[i][0]] += in[(int)double_A[i][1]] * double_A[i][2];
	}
}
void c_av(int n, s_complex* in, s_complex* out)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		out[i] = s_Complex(0, 0);
	}
	for (i = 0; i < tlen; i++)
	{
		out[(int)creal(complex_A[i][0])] += in[(int)creal(complex_A[i][1])] * complex_A[i][2];
	}
}
void dsaupd(int n, int nev, double* Evals, double** Evecs)
{
	int ido = 0;
	char bmat[2] = "I";
	char which[3] = "SM";
	double tol = 0.0;
	double* resid;
	resid = new double[n];
	int ncv = 4 * nev;
	if (ncv > n)
	{
		ncv = n;
	}
	double* v;
	int ldv = n;
	v = new double[ldv * ncv];
	int* iparam;
	iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = 3 * n;
	iparam[6] = 1;
	int* ipntr;
	ipntr = new int[11];
	double* workd;
	workd = new double[3 * n];
	double* workl;
	workl = new double[ncv * (ncv + 8)];
	int lworkl = ncv * (ncv + 8);
	int info = 0;
	int rvec = 1;
	int* select;
	select = new int[ncv];
	double* d;
	d = new double[2 * ncv];
	double sigma;
	int ierr;
	do
	{
		dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
			&ncv, v, &ldv, iparam, ipntr, workd, workl,
			&lworkl, &info);
		if ((ido == 1) || (ido == -1))
		{
			av(n, workd + ipntr[0] - 1, workd + ipntr[1] - 1);
		}
	} while ((ido == 1) || (ido == -1));
	if (info < 0)
	{
		cout << "Error with dsaupd, info = " << info << "\n";
		cout << "Check documentation in dsaupd\n\n";
	}
	else
	{
		char All[] = "All";
		dseupd_(&rvec, All, select, d, v, &ldv, &sigma, bmat,
			&n, which, &nev, &tol, resid, &ncv, v, &ldv,
			iparam, ipntr, workd, workl, &lworkl, &ierr);
		if (ierr != 0)
		{
			cout << "Error with dseupd, info = " << ierr << "\n";
			cout << "Check the documentation of dseupd.\n\n";
		}
		else if (info == 1)
		{
			cout << "Maximum number of iterations reached.\n\n";
		}
		else if (info == 3)
		{
			cout << "No shifts could be applied during implicit\n";
			cout << "Arnoldi update, try increasing NCV.\n\n";
		}
		for (int i = 0; i < nev; i++)
		{
			Evals[i] = d[i];
		}
		for (int i = 0; i < nev; i++)
		{
			for (int j = 0; j < n; j++)
			{
				Evecs[j][i] = v[i * n + j];
			}
		}
		delete resid;
		delete v;
		delete iparam;
		delete ipntr;
		delete workd;
		delete workl;
		delete select;
		delete d;
	}
}
void cnaupd(int n, int nev, s_complex* Evals, s_complex** Evecs)
{
	int ido = 0;
	char bmat[2] = "I";
	char which[3] = "SM";
	double tol = 0.0;
	s_complex* resid;
	resid = new s_complex[n];
	int ncv = 4 * nev;
	if (ncv > n)
	{
		ncv = n;
	}
	s_complex* v;
	int ldv = n;
	v = new s_complex[ldv * ncv];
	int* iparam;
	iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = 3 * n;
	iparam[6] = 1;
	int* ipntr;
	ipntr = new int[11];
	s_complex* workd;
	workd = new s_complex[3 * n];
	s_complex* workl;
	int lworkl = ncv * (3 * ncv + 5);
	workl = new s_complex[lworkl];
	double* rwork;
	rwork = new double[ncv];
	int info = 0;
	do
	{
		cnaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
			&ncv, v, &ldv, iparam, ipntr, workd, workl,
			&lworkl, rwork, &info);
		if ((ido == 1) || (ido == -1))
		{
			c_av(n, workd + ipntr[0] - 1, workd + ipntr[1] - 1);
		}
	} while ((ido == 1) || (ido == -1));
	int rvec = 0;
	int* select;
	select = new int[ncv];
	s_complex* d;
	d = new s_complex[n];
	s_complex* z;
	z = new s_complex[n * nev];
	int ldz = n;
	double sigmar, sigmai;
	int ierr;
	if (info < 0)
	{
		cout << "Error with cnaupd, info = " << info << "\n";
		cout << "Check documentation in cnaupd\n\n";
	}
	else
	{
		char All[] = "All";
		cneupd_(&rvec, All, select, d, v, &ldv, &sigmar, &sigmai,
			bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
			iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);
		if (ierr != 0)
		{
			cout << "Error with cnaupd, info = " << ierr << "\n";
			cout << "Check the documentation of cnaupd.\n\n";
		}
		else if (info == 1)
		{
			cout << "Maximum number of iterations reached.\n\n";
		}
		else if (info == 3)
		{
			cout << "No shifts could be applied during implicit\n";
			cout << "Arnoldi update, try increasing NCV.\n\n";
		}
		for (int i = 0; i < nev; i++)
		{
			Evals[i] = s_Complex(creal(d[nev - 1 - i]), cimag(d[nev - 1 - i]));
		}
		for (int i = 0; i < nev; i++)
		{
			for (int j = 0; j < n; j++)
			{
				Evecs[j][i] = s_Complex(creal(v[i * n + j]), cimag(v[i * n + j]));
			}
		}
		delete resid;
		delete v;
		delete iparam;
		delete ipntr;
		delete workd;
		delete workl;
		delete select;
		delete d;
	}
}

void dDiff(int n, int nev, double* Evals, double** Evecs, double** res)
{
	for (int i = 0; i < nev; i++)
	{
		double* Tv = new double[n];
		double* Ev = new double[n];
		for (int j = 0; j < n; j++)
		{
			Ev[j] = Evals[i] * Evecs[j][i];
			Tv[j] = 0;
			for (int k = 0; k < n; k++)
			{
				Tv[j] += Evecs[k][i] * double_A[k + j * n][2];
			}
			res[j][i] = Tv[j] - Ev[j];
		}
	}
}

void cDiff(int n, int nev, s_complex* Evals, s_complex** Evecs, s_complex** res)
{
	for (int i = 0; i < nev; i++)
	{
		s_complex* Tv = new s_complex[n];
		s_complex* Ev = new s_complex[n];
		for (int j = 0; j < n; j++)
		{
			Ev[j] = Evals[i] * Evecs[j][i];
			Tv[j] = s_Complex(0, 0);
			for (int k = 0; k < n; k++)
			{
				Tv[j] += Evecs[k][i] * complex_A[k + j * n][2];
			}
			res[j][i] = Tv[j] - Ev[j];
		}
	}

}

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "RUS");
	double t_ev;
	int n, nev;
	vector<double> vres;
	int max1 = 2100, min1 = 1000, counter = 100;
	//ввод-информация
	if (argv[1][0] == 'a') cout << "Выбран вариант а" << endl;
	else cout << "Выбран вариант b" << endl;

	cout << "Выводить матрицы (1 - yes, 0 - no)?" << endl;
	cin >> fact;

	if (fact)
	{
		cout << "Размер выводимой матрицы/вектора (3) (часть от 2-15 или 0(полностью)):" << endl;
		cin >> numbermatr;
		if (numbermatr == 0) real_size = true;
	}


	cout << "Значение N для min матрицы/вектора (1000):" << endl;
	cin >> min1;

	cout << "Значение N для max матрицы/вектора (2100):" << endl;
	cin >> max1;

	cout << "Увеличивать размер матрицы на (100) элементов:" << endl;
	cin >> counter;

	for (int n = min1; n < max1 + 1; n += counter)
	{
		nev = n - 2;
		cout << "\n№" << (((n - min1) + counter) / counter) << "  Размер N: " << "" << n << endl;

		if (real_size) numbermatr = n;  //реальный размер выводимых матриц/вектора


		if (argv[1][0] == 'a')
		{
			int i, j;
			double** Evecs, * Evals, ** nv;
			tlen = n * n;
			double_A = new double* [tlen];
			for (i = 0; i < tlen; i++)
			{
				double_A[i] = new double[3];
			}

			tlen = 0;
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					double_A[tlen][0] = i;
					double_A[tlen][1] = j;
					double_A[tlen][2] = i + j + fabs(i - j);
					tlen++;
				}
			}
			//generation
			/*if (fact) {
				cout << "\nГенерация матрицы ( " << n << "x" << n << " ):" << endl;
				for (int i = 0; i < numbermatr; i++) {
					for (int j = 0; j < numbermatr; j++) {
						cout << double_A[i + j * numbermatr][2] << "\t";
					}
					cout << endl;
				}
			}*/

			Evals = new double[nev];
			Evecs = new double* [n];
			nv = new double* [n];
			for (i = 0; i < n; i++)
			{
				Evecs[i] = new double[nev];
				nv[i] = new double[nev];
			}
			float t1 = clock();
			dsaupd(n, nev, Evals, Evecs);
			t_ev = (double)(clock() - t1) / CLOCKS_PER_SEC;
			//dDiff(n, nev, Evals, Evecs, nv);
/*
			//--------------------------------
			for (int i = 0; i < nev; i++)
			{
				double* Tv = new double[n];
				double* Ev = new double[n];
				for (int j = 0; j < n; j++)
				{
					Ev[j] = Evals[i] * Evecs[j][i];
					Tv[j] = 0;
					for (int k = 0; k < n; k++)
					{
						Tv[j] += Evecs[k][i] * double_A[k + j * n][2];
					}
					nv[j][i] = Tv[j] - Ev[j];

				}
			}

			if (fact) {

				cout << "\nМатрица A ( " << n << "x" << n << " ):" << endl;
				//printMatrix
				for (int i = 0; i < numbermatr; i++) {
					for (int j = 0; j < numbermatr; j++) {
						cout << double_A[i + j * numbermatr][2] << "\t";
					}
					cout << endl;
				}

				//printvalue and vector
				for (int i = 0; i < numbermatr - 2; i++) {
					cout << "\n'Собственное значение' " << i + 1 << ": " << Evals[i] << endl;

					cout << "'Собственный вектор' " << i + 1 << ": (";
					for (j = 0; j < numbermatr; j++)
					{
						cout << Evecs[j][i] << "  ";
					}
					cout << ")" << endl;

					cout << "Невязка (" << i + 1 << "): (";
					for (int j = 0; j < numbermatr; j++)
					{
						cout << nv[j][i] << "  ";
					}
					cout << ")" << endl;
				}
			}*/
			//------------------------------------

		}
		else
		{
			int i, j;
			s_complex** Evecs, * Evals, ** nv;
			tlen = n * n;
			complex_A = new s_complex * [tlen];
			for (i = 0; i < tlen; i++)
			{
				complex_A[i] = new s_complex[3];
			}
			tlen = 0;
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					complex_A[tlen][0] = s_Complex(i, 0);
					complex_A[tlen][1] = s_Complex(j, 0);
					complex_A[tlen][2] = s_Complex(i + j, fabs(2 * i - j));
					tlen++;
				}
			}
			//generation
			/*if (fact) {
				cout << "\nГенерация матрицы ( " << n << "x" << n << " ):" << endl;
				for (int i = 0; i < numbermatr; i++) {
					for (int j = 0; j < numbermatr; j++) {
						printf("%.2f %+.2fi\t", creal(complex_A[i + j * numbermatr][2]), cimag(complex_A[i + j * numbermatr][2]));
					}
					cout << endl;
				}
			}*/


			Evals = new s_complex[nev];
			Evecs = new s_complex * [n];
			nv = new s_complex * [n];
			for (i = 0; i < n; i++)
			{
				Evecs[i] = new s_complex[nev];
				nv[i] = new s_complex[nev];
			}
			float t1 = clock();
			cnaupd(n, nev, Evals, Evecs);
			t_ev = (double)(clock() - t1) / CLOCKS_PER_SEC;
			//cDiff(n, nev, Evals, Evecs, nv);
/*
			//----------------------------------
			for (int i = 0; i < nev; i++)
			{
				s_complex* Tv = new s_complex[n];
				s_complex* Ev = new s_complex[n];
				for (int j = 0; j < n; j++)
				{
					Ev[j] = Evals[i] * Evecs[j][i];
					Tv[j] = s_Complex(0, 0);
					for (int k = 0; k < n; k++)
					{
						Tv[j] += Evecs[k][i] * complex_A[k + j * n][2];
					}
					nv[j][i] = Tv[j] - Ev[j];
				}
			}


			if (fact) {

				cout << "\nМатрица A ( " << n << "x" << n << " ):" << endl;
				//printMatrix
				for (int i = 0; i < numbermatr; i++) {
					for (int j = 0; j < numbermatr; j++) {
						printf("%.2f %+.2fi\t", creal(complex_A[i + j * numbermatr][2]), cimag(complex_A[i + j * numbermatr][2]));
					}
					cout << endl;
				}

				//printvalue and vector with Discrepancy
				for (int i = 0; i < numbermatr - 2; i++) {
					cout << "\n'Собственное значение' " << i + 1 << " : ";
					printf("%.2f %+.2fi  ", creal(Evals[i]), cimag(Evals[i]));

					cout << "\n'Собственный вектор' " << i + 1 << ": (";
					for (j = 0; j < numbermatr; j++)
					{
						printf("%.2f %+.2fi  ", creal(Evecs[j][i]), cimag(Evecs[j][i]));
					}
					cout << ")" << endl;

					cout << "Невязка (" << i + 1 << "): (";
					for (int j = 0; j < numbermatr; j++)
					{
						printf("%.2f %+.2fi  ", creal(nv[j][i]), cimag(nv[j][i]));
					}
					cout << ")" << endl;
				}
			}
*/
			//----------------------------------

		}
		cout << "Time of EV ( " << n << "x" << n << " ): " << t_ev << " s" << endl;
		vres.push_back(t_ev);
	}
	cout << endl;

	for (int i = 0; i < vres.size(); i++)
		cout << vres[i] << endl;
	return 0;
}
