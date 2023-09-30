#include <iostream>
#include <complex.h>
#include <ctime>
#include <cstring>
#include <cmath>
#ifdef __cplusplus
extern "C"
{
#endif
#include <cblas.h>
#include <lapacke.h>
#ifdef __cplusplus
}
#endif

using namespace std;

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "RUS");
	if (argc < 2)
	{
		cout << "������������ �� ������ ����� �������: + ��. ��� ����. �������� � ���������� �� ����.'" << endl;
		return 1;
	}
	if (argc > 2)
	{
		cout << "���� ������� ����� ������ ���������" << endl;
		return 2;
	}
	if (strlen(argv[1]) > 1 &&
		(argv[1][0] == 'e') || (argv[1][0] == 'y') || (argv[1][0] == 'u') || (argv[1][0] == 'i') || (argv[1][0] == 'o') || (argv[1][0] == 'a') || (argv[1][0] == 'q') || (argv[1][0] == 'j'))
	{
		true;
	}
	else {
		cout << "������������ �������� �����" << endl;
		return 3;
	}
	int n = 3;

	if ((argv[1][0] == 'e') || (argv[1][0] == 'y') || (argv[1][0] == 'u') || (argv[1][0] == 'i') || (argv[1][0] == 'o') || (argv[1][0] == 'a') || (argv[1][0] == 'q') || (argv[1][0] == 'j'))
	{
		cout << "������ ����� ������� ������� " << argv[1][0] << endl;
		float* mA;
		float* mB;
		float* v;
		float* av;
		float* mv;
		float* mm;
		float* mT;
		cout << "������� ������ ������� (�������� ��. ���� �� 1 �� 100): " << endl;
		cin >> n;
		mA = new float[n * n];
		mB = new float[n * n];
		v = new float[n];
		av = new float[n];
		mv = new float[n];
		mm = new float[n * n];
		mT = new float[n * n];
		srand(time(0));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				mA[i + j * n] = rand() % 100 - 0 + 1;
				mB[i + j * n] = rand() % 100 - 0 + 1;
			}
			v[i] = rand() % 100 - 0 + 1;
		}
		cout << "������� A: " << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mA[i + j * n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "������� B: " << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mB[i + j * n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "������ ��������� ������: " << endl;
		for (int i = 0; i < n; i++)
		{
			cout << v[i] << " ";
		}
		cout << endl << endl;
		cblas_sgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, mB, n, v, 1, 0.0, mv, 1);  			// B * ������ ��������� ������
		cblas_sscal(n, 10.0, v, 1);									//10 * ������ ��������� ������
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mA, n, mB, n, 0.0, mm, n); // A * B
		cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, n, n, n, 1.0, mA, n, mB, n, 0.0, mT, n); // A * B
		int matrix_layout = CblasColMajor;
		char jobvl = 'V';
/*(�������) ������*1
 = 'N': ����� ����������� ������� A �� �����������;
 = 'V': ����������� ����� ����������� ������� A.*/
		char jobvr = 'V';
/*JOBVR (�������)*1
 = 'N': ������ ����������� ������� A �� �����������;
 = 'V': ����������� ������ ����������� ������� A.*/
		lapack_int N = n;
		float* a = mA;
		lapack_int lda = n;
/*(�������) ����� �����
 ������� ������ ������� A. LDA >=
 ����.(1,N).*/
		float* wr = new float[n]; 
//������������ ������, ����������� (N)
		float* wi = new float[n]; 
/*������������ ������, ����������� (N) WR � WI �������� �������������� � ������ ����� ����������� ����������� �������� ��������������. ���������� �����������
*���� ����������� �������� ���������� ���������������, ������
����������� �������� ������� ����� ������������� ������ �����.*/
		float* vl = new float[n * n]; 
/*������������ ������, ����������� (LDVL,N)
 ���� JOBVL = 'V', ����� ����������� ������� u(j)
����������� ���� �� ������ � �������� VL �
��� �� �������, ��� � �� ����������� ��������. ���� JOBVL =
 'N', �� VL ������ ���. ���� j-� ����������� ��������
������������, �� u(j) = VL(:,j), j-� ������� VL.
 ���� j-� � (j+1)-� ����������� �������� �������� ����������
����������� ����, �� u(j) = VL(:,j) + i*VL(:,j+1)
�
u(j+1) = VL(:,j) = i*VL(:,j+1).*/
		lapack_int ldvl = n;
/*����� �����
 ������� ������ ������� VL. LDVL >= 1;
���� JOBVL = 'V', �� LDVL >= N.*/
		float* vr = new float[n * n];
/*
 ������� ��������� ������� VR. LDVR >= 1;
���� JOBVR = 'V', �� LDVR >= N.*/
		lapack_int ldvr = n;
		lapack_int inf;
		inf = LAPACKE_sgeev(matrix_layout, jobvl, jobvr, N, a, lda, wr, wi, vl, ldvl, vr, ldvr);
/*SGEEV - ��������� ��� N-��-N ������������ �������������� ������� A
����������� �������� �, �������������, ����� �/��� ������
����������� �������*/
		cout << "10 * ������ ��������� ������: " << endl;
		for (int i = 0; i < n; i++)
		{
			cout << v[i] << " ";
		}
		cout << endl << endl;
		cout << "������� B * ������ ��������� ������: " << endl;
		for (int i = 0; i < n; i++)
		{
			cout << mv[i] << " ";
		}
		cout << endl << endl;
		cout << "������� A * ������� B: " << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mm[i + j * n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "������� T(A) * T(B): " << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mT[i + j * n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
		int i = 0;
		int j;
		while (i < n)
		{
			if (wi[i] == 0.0)
			{
				cout << "WR - ������. �. ���. ������. ����. " << i << ": " << wr[i] << endl;
				cout << "����� ������. ������� u(j) VL " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vl[j + i * n];
				}
				cout << endl;
				cout << "������� ��������� ������� VR " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vr[j + i * n];
				}
				cout << endl << endl;
				i++;
			}
			else
			{
				cout << "������/������. �. " << i << ": " << wr[i] << (wi[i] > 0 ? "+" : "") << wi[i] << "i" << endl;
				cout << "����� ������. ������� u(j) VL " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vl[j + i * n] << (vl[j + (i + 1) * n] > 0 ? "+" : "") << vl[j + (i + 1) * n] << "i";
				}
				cout << endl;
				cout << "������� ��������� ������� VR " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vr[j + i * n] << (vr[j + (i + 1) * n] > 0 ? "+" : "") << vr[j + (i + 1) * n] << "i";
				}
				cout << endl << endl;
				cout << "������/������. �. " << i << ": " << wr[i] << (wi[i] > 0 ? "-" : "+") << fabs(wi[i]) << "i" << endl;
				cout << "����� ������. ������� u(j) VL " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vl[j + i * n] << (vl[j + (i + 1) * n] > 0 ? "-" : "+") << fabs(vl[j + (i + 1) * n]) << "i";
				}
				cout << endl;
				cout << "������� ��������� ������� VR " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vr[j + i * n] << (vr[j + (i + 1) * n] > 0 ? "-" : "+") << fabs(vr[j + (i + 1) * n]) << "i";
				}
				cout << endl << endl;
				i += 2;
			}
		}
		delete[] mA;
		delete[] mB;
		delete[] v;
		delete[] mv;
		delete[] mm;
		delete[] a;
		delete[] wr;
		delete[] wi;
		delete[] vl;
		delete[] vr;
	}
	/* ����������� ����� - "������" ������ ������� - ������ ������
	else
	{
		float _Complex* mA;
		float _Complex* mB;
		float _Complex* v;
		float _Complex* av;
		float _Complex* mv;
		float _Complex* mm;
		mA = new float _Complex[n * n];
		mB = new float _Complex[n * n];
		v = new float _Complex[n];
		av = new float _Complex[n];
		mv = new float _Complex[n];
		mm = new float _Complex[n * n];
		srand(time(0));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				mA[i + j * n] = (float)(rand() % 100 - 50) + (float)(rand() % 100 - 50) * _Complex_I;
				mB[i + j * n] = (float)(rand() % 100 - 50) + (float)(rand() % 100 - 50) * _Complex_I;
			}
			v[i] = (float)(rand() % 100 - 50) + (float)(rand() % 100 - 50) * _Complex_I;
		}
		cout << "Matrix A:" << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << creal(mA[i + j * n]) << (cimag(mA[i + j * n]) > 0 ? "+" : "") << cimag(mA[i + j * n]) << "i ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Matrix B:" << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << creal(mB[i + j * n]) << (cimag(mB[i + j * n]) > 0 ? "+" : "") << cimag(mB[i + j * n]) << "i ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Vector:" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << creal(v[i]) << (cimag(v[i]) > 0 ? "+" : "") << cimag(v[i]) << "i ";
		}
		cout << endl << endl;
		float _Complex alpha = (float)(1) + (float)(0) * _Complex_I;
		float _Complex beta = (float)(0) + (float)(0) * _Complex_I;
		cblas_cgemv(CblasColMajor, CblasNoTrans, n, n, &alpha, mA, n, v, 1, &beta, mv, 1);
		float _Complex constKoef = (float)(15) + (float)(-6) * _Complex_I;
		cblas_cscal(n, &constKoef, v, 1);
		cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha, mA, n, mB, n, &beta, mm, n);
		int matrix_layout = CblasColMajor;
		char jobvl = 'V';
		char jobvr = 'V';
		lapack_int N = n;
		lapack_complex_float* a = mA;
		lapack_int lda = n;
		lapack_complex_float* w = new lapack_complex_float[n];
		lapack_complex_float* vl = new lapack_complex_float[n * n];
		lapack_int ldvl = n;
		lapack_complex_float* vr = new lapack_complex_float[n * n];
		lapack_int ldvr = n;
		lapack_int inf;
		inf = LAPACKE_cgeev(matrix_layout, jobvl, jobvr, N, a, lda, w, vl, ldvl, vr, ldvr);
		cout << "(15-6i) * Vector:" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << creal(v[i]) << (cimag(v[i]) > 0 ? "+" : "") << cimag(v[i]) << "i ";
		}
		cout << endl << endl;
		cout << "Matrix A * Vector:" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << creal(mv[i]) << (cimag(mv[i]) > 0 ? "+" : "") << cimag(mv[i]) << "i ";
		}
		cout << endl << endl;
		cout << "Matrix A * Matrix B:" << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << creal(mm[i + j * n]) << (cimag(mm[i + j * n]) > 0 ? "+" : "") << cimag(mm[i + j * n]) << "i ";
			}
			cout << endl;
		}
		cout << endl;
		int i = 0;
		int j;
		while (i < n)
		{
			cout << "EV " << i << ": " << creal(w[i]) << (cimag(w[i]) > 0 ? "+" : "") << cimag(w[i]) << "i" << endl;
			cout << "VL " << i << ":";
			for (j = 0; j < n; j++)
			{
				cout << " " << creal(vl[i + j * n]) << (cimag(vl[i + j * n]) > 0 ? "+" : "") << cimag(vl[i + j * n]) << "i";
			}
			cout << endl;
			cout << "VR " << i << ":";
			for (j = 0; j < n; j++)
			{
				cout << " " << creal(vr[i + j * n]) << (cimag(vr[i + j * n]) > 0 ? "+" : "") << cimag(vr[i + j * n]) << "i";
			}
			cout << endl << endl;
			i++;
		}
		delete[] mA;
		delete[] mB;
		delete[] v;
		delete[] mv;
		delete[] mm;
		delete[] a;
		delete[] w;
		delete[] vl;
		delete[] vr;

	}*/

	return 0;
}
