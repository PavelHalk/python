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
		cout << "Определяется по второй букве фамилии: + гл. или согл. Напишите её мальенькой по англ.'" << endl;
		return 1;
	}
	if (argc > 2)
	{
		cout << "Было введено более одного параметра" << endl;
		return 2;
	}
	if (strlen(argv[1]) > 1 &&
		(argv[1][0] == 'e') || (argv[1][0] == 'y') || (argv[1][0] == 'u') || (argv[1][0] == 'i') || (argv[1][0] == 'o') || (argv[1][0] == 'a') || (argv[1][0] == 'q') || (argv[1][0] == 'j'))
	{
		true;
	}
	else {
		cout << "Недопустимый параметр введён" << endl;
		return 3;
	}
	int n = 3;

	if ((argv[1][0] == 'e') || (argv[1][0] == 'y') || (argv[1][0] == 'u') || (argv[1][0] == 'i') || (argv[1][0] == 'o') || (argv[1][0] == 'a') || (argv[1][0] == 'q') || (argv[1][0] == 'j'))
	{
		cout << "Вторая буква фамилии гласная " << argv[1][0] << endl;
		float* mA;
		float* mB;
		float* v;
		float* av;
		float* mv;
		float* mm;
		float* mT;
		cout << "Впишите размер матрицы (значения эл. ранд от 1 до 100): " << endl;
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
		cout << "Матрица A: " << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mA[i + j * n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Матрица B: " << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mB[i + j * n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Вектор свободных членов: " << endl;
		for (int i = 0; i < n; i++)
		{
			cout << v[i] << " ";
		}
		cout << endl << endl;
		cblas_sgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, mB, n, v, 1, 0.0, mv, 1);  			// B * Вектор свободных членов
		cblas_sscal(n, 10.0, v, 1);									//10 * Вектор свободных членов
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, mA, n, mB, n, 0.0, mm, n); // A * B
		cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, n, n, n, 1.0, mA, n, mB, n, 0.0, mT, n); // A * B
		int matrix_layout = CblasColMajor;
		char jobvl = 'V';
/*(входной) СИМВОЛ*1
 = 'N': левые собственные векторы A не вычисляются;
 = 'V': вычисляются левые собственные векторы A.*/
		char jobvr = 'V';
/*JOBVR (входной)*1
 = 'N': правые собственные векторы A не вычисляются;
 = 'V': вычисляются правые собственные векторы A.*/
		lapack_int N = n;
		float* a = mA;
		lapack_int lda = n;
/*(входное) ЦЕЛОЕ число
 Ведущий размер массива A. LDA >=
 макс.(1,N).*/
		float* wr = new float[n]; 
//ВЕЩЕСТВЕННЫЙ массив, размерность (N)
		float* wi = new float[n]; 
/*ВЕЩЕСТВЕННЫЙ массив, размерность (N) WR и WI содержат действительную и мнимую части вычисленных собственных значений соответственно. Комплексно сопряженные
*пары собственных значений появляются последовательно, причем
собственное значение сначала имеет положительную мнимую часть.*/
		float* vl = new float[n * n]; 
/*ВЕЩЕСТВЕННЫЙ массив, размерность (LDVL,N)
 Если JOBVL = 'V', левые собственные векторы u(j)
сохраняются один за другим в столбцах VL в
том же порядке, что и их собственные значения. Если JOBVL =
 'N', на VL ссылки нет. Если j-е собственное значение
вещественное, то u(j) = VL(:,j), j-й столбец VL.
 Если j-е и (j+1)-е собственные значения образуют комплексно
сопряженную пару, то u(j) = VL(:,j) + i*VL(:,j+1)
и
u(j+1) = VL(:,j) = i*VL(:,j+1).*/
		lapack_int ldvl = n;
/*целое число
 Ведущий размер массива VL. LDVL >= 1;
если JOBVL = 'V', то LDVL >= N.*/
		float* vr = new float[n * n];
/*
 Ведущее измерение массива VR. LDVR >= 1;
если JOBVR = 'V', то LDVR >= N.*/
		lapack_int ldvr = n;
		lapack_int inf;
		inf = LAPACKE_sgeev(matrix_layout, jobvl, jobvr, N, a, lda, wr, wi, vl, ldvl, vr, ldvr);
/*SGEEV - вычислите для N-на-N вещественной несимметричной матрицы A
собственные значения и, необязательно, левый и/или правый
собственные векторы*/
		cout << "10 * Вектор свободных членов: " << endl;
		for (int i = 0; i < n; i++)
		{
			cout << v[i] << " ";
		}
		cout << endl << endl;
		cout << "Матрица B * Вектор свободных членов: " << endl;
		for (int i = 0; i < n; i++)
		{
			cout << mv[i] << " ";
		}
		cout << endl << endl;
		cout << "Матрица A * Матрица B: " << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mm[i + j * n] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "Матрица T(A) * T(B): " << endl;
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
				cout << "WR - действ. ч. выч. собств. знач. " << i << ": " << wr[i] << endl;
				cout << "левые собств. векторы u(j) VL " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vl[j + i * n];
				}
				cout << endl;
				cout << "Ведущее измерение массива VR " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vr[j + i * n];
				}
				cout << endl << endl;
				i++;
			}
			else
			{
				cout << "мнимая/действ. ч. " << i << ": " << wr[i] << (wi[i] > 0 ? "+" : "") << wi[i] << "i" << endl;
				cout << "левые собств. векторы u(j) VL " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vl[j + i * n] << (vl[j + (i + 1) * n] > 0 ? "+" : "") << vl[j + (i + 1) * n] << "i";
				}
				cout << endl;
				cout << "Ведущее измерение массива VR " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vr[j + i * n] << (vr[j + (i + 1) * n] > 0 ? "+" : "") << vr[j + (i + 1) * n] << "i";
				}
				cout << endl << endl;
				cout << "мнимая/действ. ч. " << i << ": " << wr[i] << (wi[i] > 0 ? "-" : "+") << fabs(wi[i]) << "i" << endl;
				cout << "левые собств. векторы u(j) VL " << i << ":";
				for (j = 0; j < n; j++)
				{
					cout << " " << vl[j + i * n] << (vl[j + (i + 1) * n] > 0 ? "-" : "+") << fabs(vl[j + (i + 1) * n]) << "i";
				}
				cout << endl;
				cout << "Ведущее измерение массива VR " << i << ":";
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
	/* Комплексные числа - "Чурсин" вторая гласная - значит лишнее
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
