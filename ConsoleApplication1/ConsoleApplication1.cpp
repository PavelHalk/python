//lu decomposition C++

#include<iostream>
#include<iomanip>

//size of matrix
# define N 4
typedef float matrix[N][N];
matrix l, u, a;

using namespace std;

float b[N], x[N], v[N];

//upper triangular matrix

void urow(int i)
{
	float s;
	int j, k;
	for (j = i; j < N; j++)
	{
		s = 0;
		for (k = 0; k < i; k++)
			s += u[k][j] * l[i][k];
		u[i][j] = a[i][j] - s;
	}
}

//lower triangular matrix

void lcol(int j)
{
	float s;
	int i, k;
	for (i = j + 1; i < N; i++)
	{
		s = 0;
		for (k = 0; k <= i - 1; k++)
			s += u[k][j] * l[i][k];
		l[i][j] = (a[i][j] - s) / u[j][j];
	}
}
//function to display the matrix

void printmat(matrix x)
{
	int i, j;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			cout << setw(8) << setprecision(4) << x[i][j];
		cout << endl;
	}
}
int main()
{
	int i, j, m;
	float s;
	cout << "Enter the elements of augmented matrix row wise" << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			cin >> a[i][j];
		cin >> b[i];
	}
	cout << "\nInput matrix: " << endl;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
			cout << a[i][j] << '	';
		cout << b[i] << endl;
	}
	cout << endl;
	for (i = 0; i < N; i++)
		l[i][i] = 1.0;
	for (m = 0; m < N; m++)
	{
		urow(m);
		if (m < (N - 1))
			lcol(m);
	}
	cout << setw(14) << "U \n" << endl;
	printmat(u);
	cout << endl;
	cout << setw(14) << "L \n" << endl;
	printmat(l);

	for (i = 0; i < N; i++)
	{
		s = 0;
		for (j = 0; j <= i - 1; j++)
			s += l[i][j] * v[j];
		v[j] = b[i] - s;
	}
	for (i = N - 1; i >= 0; i--)
	{
		s = 0;
		for (j = i + 1; j < N; j++)
			s += u[i][j] * x[j];
		x[i] = (v[i] - s) / u[i][i];
	}
	cout << "The solution is: " << endl;
	for (i = 0; i < N; i++)
		cout << "x[" << setw(3) << i + 1 << "] =" << setw(6) << setprecision(4) << x[i] << endl;
	return 0;
}