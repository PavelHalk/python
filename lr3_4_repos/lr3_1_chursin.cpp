#include <iostream>
#include <cstring>
#include <complex.h>
#include <ctime>
#include <vector>

using namespace std;

bool fact = false, real_size = false; //print matrix
int numbermatr = 3; //выводимая часть

void LU_real(int n, double* a, double* l, double* u)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                l[i + j * n] = 1;
            }
            else
            {
                l[i + j * n] = 0;
            }
            u[i + j * n] = 0;
        }
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double tmp = 0;
            if (i <= j)
            {
                for (int k = 0; k < i; k++)
                {
                    tmp += l[i + k * n] * u[k + j * n];
                }
                u[i + j * n] = a[i + j * n] - tmp;
            }
            else
            {
                for (int k = 0; k < j; k++)
                {
                    tmp += l[i + k * n] * u[k + j * n];
                }
                l[i + j * n] = (a[i + j * n] - tmp) / u[j + j * n];
            }
        }
    }
}
void LU_Solver_real(int n, double* a, double* b)
{
    double* l = new double[n * n];
    double* u = new double[n * n];
    double* result = new double[n * n];
    LU_real(n, a, l, u);
    double* t = new double[n];
    for (int i = 0; i < n; i++)
    {
        double tmp = 0;
        t[i] = 0;
        for (int j = 0; j < i; j++)
        {
            tmp += t[j] * l[i + j * n];
        }
        t[i] = (b[i] - tmp) / l[i + i * n];
    }
    for (int i = n - 1; i >= 0; i--)
    {
        double tmp = 0;
        b[i] = 0;
        for (int j = n - 1; j >= 0; j--)
        {
            tmp += b[j] * u[i + j * n];
        }
        b[i] = (t[i] - tmp) / u[i + i * n];
    }
    
    if (fact) {
        cout << "\nЧасть матрицы l (" << numbermatr << "):" << endl;
        /*printMatrix*/
        for (int i = 0; i < numbermatr; i++) {
            for (int j = 0; j < numbermatr; j++) {
                if (i != 0 && !real_size) cout << l[i * numbermatr + j - 1] << "\t";
                else cout << l[i * numbermatr + j] << "\t";
            }
            cout << endl;
        }

        cout << "\nЧасть матрицы u (" << numbermatr << "):" << endl;
        /*printMatrix*/
        for (int i = 0; i < numbermatr; i++) {
            for (int j = 0; j < numbermatr; j++) {
                if (i != 0 && !real_size) cout << u[i * numbermatr + j - 1] << "\t";
                else cout << u[i * numbermatr + j] << "\t";
            }
            cout << endl;
        }
    }
    
    /*cout << "\nЧасть матрицы l*u (" << numbermatr << "):" << endl;
        /*printMatrix      
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i * j] = 0;
            for (int k = 0; k < n; k++) {
                result[i*j] += ( u[i * k]  * l[k * j]);
            }
        }
    }
    
        /*printMatrix
        for (int i = 0; i < numbermatr; i++) {
            for (int j = 0; j < numbermatr; j++) {
                if (i != 0 && !real_size) cout << result[i * numbermatr + j - 1] << "\t";
                else cout << result[i * numbermatr + j] << "\t";
            }
            cout << endl;
        }*/
    
    delete[] l;
    delete[] u;
}
void LU_complex(int n, double _Complex* a, double _Complex* l, double _Complex* u)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                l[i + j * n] = (double)1 + (double)0 * _Complex_I;
            }
            else
            {
                l[i + j * n] = (double)0 + (double)0 * _Complex_I;
            }
            u[i + j * n] = (double)0 + (double)0 * _Complex_I;
        }
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double _Complex tmp = 0;
            if (i <= j)
            {
                for (int k = 0; k < i; k++)
                {
                    tmp += l[i + k * n] * u[k + j * n];
                }
                u[i + j * n] = a[i + j * n] - tmp;
            }
            else
            {
                for (int k = 0; k < j; k++)
                {
                    tmp += l[i + k * n] * u[k + j * n];
                }
                l[i + j * n] = (a[i + j * n] - tmp) / u[j + j * n];
            }
        }
    }
}
void LU_Solver_complex(int n, double _Complex* a, double _Complex* b)
{
    double _Complex* l = new double _Complex[n * n];
    double _Complex* u = new double _Complex[n * n];
    LU_complex(n, a, l, u);
    double _Complex* t = new double _Complex[n];
    for (int i = 0; i < n; i++)
    {
        double _Complex tmp = (double)0 + (double)0 * _Complex_I;
        t[i] = (double)0 + (double)0 * _Complex_I;
        for (int j = 0; j < i; j++)
        {
            tmp += t[j] * l[i + j * n];
        }
        t[i] = (b[i] - tmp) / l[i + i * n];
    }
    for (int i = n - 1; i >= 0; i--)
    {
        double _Complex tmp = (double)0 + (double)0 * _Complex_I;
        b[i] = (double)0 + (double)0 * _Complex_I;
        for (int j = n - 1; j >= 0; j--)
        {
            tmp += b[j] * u[i + j * n];
        }
        b[i] = (t[i] - tmp) / u[i + i * n];
    }


    if (fact) {
        cout << "\nЧасть матрицы l (" << numbermatr << "):" << endl;
        /*printMatrix*/
        for (int i = 0; i < numbermatr; i++) {
            for (int j = 0; j < numbermatr; j++) {
                if (i != 0 && !real_size) printf("%.2f %+.2fi\t", creal(l[i * numbermatr + j-1]), cimag(l[i * numbermatr + j-1]));
                else printf("%.2f %+.2fi\t", creal(l[i * numbermatr + j]), cimag(l[i * numbermatr + j]));
            }
            cout << endl;
        }

        cout << "\nЧасть матрицы u (" << numbermatr << "):" << endl;
        /*printMatrix*/
        for (int i = 0; i < numbermatr; i++) {
            for (int j = 0; j < numbermatr; j++) {
                if (i != 0 && !real_size) printf("%.2f %+.2fi\t", creal(u[i * numbermatr + j-1]), cimag(u[i * numbermatr + j-1]));
                else printf("%.2f %+.2fi\t", creal(u[i * numbermatr + j]), cimag(u[i * numbermatr + j]));
            }
            cout << endl;
        }
    }


    delete[] l;
    delete[] u;
}
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
            double* a = new double[n * n];
            double* b = new double[n];
            cout << "Начало решения (LU-факторизация) dgetrf №" << (((n - min1) + counter) / counter) << endl;
            srand(time(0));
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a[i + j * n] = rand() % 10 - 5;
                }
                b[i] = rand() % 10 - 5;
            }

            if (fact) {
                cout << "\nЧасть матрицы А (" << numbermatr << "):" << endl;
                /*printMatrix*/
                for (int i = 0; i < numbermatr; i++) {
                    for (int j = 0; j < numbermatr; j++) {
                        if (i != 0 && !real_size) cout << a[i * numbermatr + j - 1] << "\t";
                        else cout << a[i * numbermatr + j] << "\t";
                    }
                    cout << endl;
                }

                cout << "\nЧасть вектора B (" << numbermatr << "):" << endl;
                /*printMatrix*/
                for (int i = 0; i < numbermatr; i++) {
                    cout << b[i] << "\t";
                }
                cout << endl;
            }

            // cout << "Формирование Random матриц завершено" << endl;
            int t1 = clock();
            LU_Solver_real(n, a, b); //решение вещественных
            v1 = (double)(clock() - t1) / CLOCKS_PER_SEC;
            cout << "\nВремя решения №" << (((n - min1) + counter) / counter) << " = " << v1 << " сек." << endl;
            delete[] a;
            delete[] b;


        }
        else
        {
            double _Complex* a = new double _Complex[n * n];
            double _Complex* b = new double _Complex[n];
            cout << "Начало решения (LU-факторизация) zgetrf №" << (((n - min1) + counter) / counter) << endl;
            srand(time(0));
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a[i + j * n] = (double)(rand() % 10 - 5) + (double)(rand() % 10 - 5) * _Complex_I;
                }
                b[i] = (double)(rand() % 10 - 5) + (float)(rand() % 10 - 5) * _Complex_I;
            }

            if (fact) {
                cout << "\nЧасть матрицы А (" << numbermatr << "):" << endl;
                /*printMatrix*/
                for (int i = 0; i < numbermatr; i++) {
                    for (int j = 0; j < numbermatr; j++) {
                        if (i != 0 && !real_size) printf("%.2f %+.2fi\t", creal(a[i * numbermatr + j-1]), cimag(a[i * numbermatr + j-1]));
                        else printf("%.2f %+.2fi\t", creal(a[i * numbermatr + j]), cimag(a[i * numbermatr + j]));
                    }
                    cout << endl;
                }

                cout << "\nЧасть вектора B (" << numbermatr << "):" << endl;
                /*printMatrix*/
                for (int i = 0; i < numbermatr; i++) {
                    printf("%.2f %+.2fi\t", creal(b[i]), cimag(b[i]));
                }
                cout << endl;
            }

            //cout << "Формирование Random матриц завершено" << endl;
            int t1 = clock();
            LU_Solver_complex(n, a, b);  //решение комплекс
            v1 = (double)(clock() - t1) / CLOCKS_PER_SEC;
            cout << "\nВремя решения №" << (((n - min1) + counter) / counter) << " = " << v1 << " сек." << endl;
            delete[] a;
            delete[] b;
        }
        vres.push_back(v1);
    }
    cout << endl;
    for (int i = 0; i < vres.size(); i++)
        cout << vres[i] << endl;
    return 0;
}
