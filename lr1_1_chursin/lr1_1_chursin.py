import numpy as np
from matrix_np_api import *
def decompose_to_LU(a):
    """Разложите матрицу коэффициентов на матрицы L и U.
    Треугольные матрицы L и U будут представлены в виде одной матрицы nxn.
    :param a: числовая матрица коэффициентов
    :return: numpy матрица LU
    """
    lu_matrix = np.matrix(np.zeros([a.shape[0], a.shape[1]]))
    n = a.shape[0]

    for k in range(n):
        # вычислить все остаточные элементы k-строки
        for j in range(k, n):
            lu_matrix[k, j] = a[k, j] - lu_matrix[k, :k] * lu_matrix[:k, j]
        # вычислить все остаточные элементы k-столбца
        for i in range(k + 1, n):
            lu_matrix[i, k] = (a[i, k] - lu_matrix[i, : k] * lu_matrix[: k, k]) / lu_matrix[k, k]

    return lu_matrix

def get_L(m):
    """Получаем треугольную L-матрицу из одной LU-матрицы
    :param m: numpy LU-матрица
    :return: numpy треугольная L-матрица
    """
    L = m.copy()
    for i in range(L.shape[0]):
            L[i, i] = 1
            L[i, i+1 :] = 0
    return np.matrix(L)


def get_U(m):
    """Получаем треугольную U-матрицу из одной LU-матрицы
    :param m: numpy LU-матрица
    :return: простая треугольная U-образная матрица
    """
    U = m.copy()
    for i in range(1, U.shape[0]):
        U[i, :i] = 0
    return U

#Напишем функцию для решения системы уравнений, исходя из полученного LU-разложения
def solve_LU(lu_matrix, b):
    """Решите систему уравнений из заданной LU-матрицы и вектора b абсолютных членов.
    :param lu_matrix: numpy LU-матрица
    :param b: числовая матрица абсолютных членов [n x 1]
    :return: числовая матрица ответов [n x 1]
    """
    # получаем опорный вектор y
    y = np.matrix(np.zeros([lu_matrix.shape[0], 1]))
    for i in range(y.shape[0]):
        y[i, 0] = b[i, 0] - lu_matrix[i, :i] * y[:i]

    # получаем вектор ответов x
    x = np.matrix(np.zeros([lu_matrix.shape[0], 1]))
    for i in range(1, x.shape[0] + 1):
        x[-i, 0] = (y[-i] - lu_matrix[-i, -i:] * x[-i:, 0] )/ lu_matrix[-i, -i]

    return x

#Для проверки правильности разложения, перемножим матрицы L и U и сравним с исходной матрицей коэффициентов
a, b = read_matrix_a_b('lu_matrix.txt')
display_matrix(a)

LU = decompose_to_LU(a)
L = get_L(LU)
U = get_U(LU)

display_matrix(L * U)