#include <iostream>
#include <complex.h>
#include <ctime>
#include <vector>
#include <cstring>
#ifdef __cplusplus
extern "C"
{
#endif //__cplusplus
//#include <cblas.h>
#include <lapack.h>
#ifdef __cplusplus
}
#endif //__cplusplus
using namespace std;

bool fact = false, real_size = false; //print matrix
int numbermatr = 3; //выводимая часть

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "RUS");

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
		double v1;
		cout << "№" << (((n - min1) + counter) / counter) << "  Размер N: " << "" << n << endl;

		if (real_size) numbermatr = n;  //реальный размер выводимых матриц/вектора

		if (argv[1][0] == 'a')
		{
			double* matr;
			matr = new double[n * n];
			int* v;
			v = new int[n];
			double* b;
			b = new double[n];
			cout << "Начало решения  dgetrf №" << (((n - min1) + counter) / counter) << endl;
			srand(time(0));
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					matr[i + j * n] = rand() % 10 - 5;
				}
				b[i] = rand() % 10 - 5;
			}
			// cout << "Формирование Random матриц завершено" << endl;
			int t1 = clock();
			char trans = 'N';
			int m = 1;
			int info;			
			// LU factor
			LAPACK_dgetrf(&n, &n, matr, &n, v, &info);
			//solve a complex system of linear equations
			LAPACK_dgetrs(&trans, &n, &m, matr, &n, v, b, &n, &info);
			v1 = (double)(clock() - t1) / CLOCKS_PER_SEC;

			if (fact) {
				cout << "\nЧасть матрицы (" << numbermatr << "):" << endl;
				/*printMatrix*/
				for (int i = 0; i < numbermatr; i++) {
					for (int j = 0; j < numbermatr; j++) {
						if (i != 0 && !real_size) cout << matr[i * numbermatr + j - 1] << "\t";
						else cout << matr[i * numbermatr + j] << "\t";
					}
					cout << endl;
				}

				cout << "\nЧасть вектора V (" << numbermatr << "):" << endl;
				/*printMatrix*/
				for (int i = 0; i < numbermatr; i++) {
					cout << v[i] << "\t";
				}
				cout << endl;

				cout << "\nЧасть вектора B (" << numbermatr << "):" << endl;
				/*printMatrix*/
				for (int i = 0; i < numbermatr; i++) {
					cout << b[i] << "\t";
				}
				cout << endl;
			}
			cout << "\nВремя решения №" << (((n - min1) + counter) / counter) << " = " << v1 << " сек." << endl;
			delete[] matr;
			delete[] v;
			delete[] b;
		}
		else
		{
			double _Complex* matr;
			matr = new double _Complex[n * n];
			int* v;
			v = new int[n];
			double _Complex* b;
			b = new double _Complex[n];
			cout << "Начало решения  zgetrf №" << (((n - min1) + counter) / counter) << endl;
			srand(time(0));
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					matr[i + j * n] = (double)(rand() % 10 - 5) + (double)(rand() % 10 - 5) * _Complex_I;
				}
				b[i] = (double)(rand() % 10 - 5) + (float)(rand() % 10 - 5) * _Complex_I;
			}
			//cout << "Формирование Random матриц завершено" << endl;
			int t1 = clock();
			int info;
			char trans = 'N';
			int m = 1;
			LAPACK_zgetrf(&n, &n, matr, &n, v, &info);
			LAPACK_zgetrs(&trans, &n, &m, matr, &n, v, b, &n, &info);
			v1 = (double)(clock() - t1) / CLOCKS_PER_SEC;

			if (fact) {
				cout << "\nЧасть матрицы (" << numbermatr << "):" << endl;
				/*printMatrix*/
				for (int i = 0; i < numbermatr; i++) {
					for (int j = 0; j < numbermatr; j++) {
						if (i != 0 && !real_size) printf("%.2f %+.2fi\t", creal(matr[i * numbermatr + j - 1]), cimag(matr[i * numbermatr + j - 1]));
						else printf("%.2f %+.2fi\t", creal(matr[i * numbermatr + j]), cimag(matr[i * numbermatr + j]));
					}
					cout << endl;
				}

				cout << "\nЧасть вектора v (" << numbermatr << "):" << endl;
				/*printMatrix*/
				for (int i = 0; i < numbermatr; i++) {
					printf("%.2f %+.2fi\t", creal(v[i]), cimag(v[i]));
				}
				cout << endl;

				cout << "\nЧасть вектора B (" << numbermatr << "):" << endl;
				/*printMatrix*/
				for (int i = 0; i < numbermatr; i++) {
					printf("%.2f %+.2fi\t", creal(b[i]), cimag(b[i]));
				}
				cout << endl;
			}
			cout << "\nВремя решения №" << (((n - min1) + counter) / counter) << " = " << v1 << " сек." << endl;

			delete[] matr;
			delete[] v;
			delete[] b;
		}
		vres.push_back(v1);
	}
	for (int i = 0; i < vres.size(); i++)
		cout << vres[i] << endl;
	return 0;
}
