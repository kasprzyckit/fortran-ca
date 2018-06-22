import sys
import numpy as np
from src.matrix_lib import matrix_lib as ml

def test_mm(N):
    first = np.ndarray(shape=(N, N), dtype=float, order='F')
    second = np.ndarray(shape=(N, N), dtype=float, order='F')

    for i in range(N):
        for j in range(N):
            first[i][j] = 2
            second[i][j] = 3

    multiply = ml.mm(first, second)
    print(multiply)

def test_gaussian(N):
    h2 = 1.0 / (N * N)
    P1 = 1.0 / h2
    P2 = -2.0 / h2

    A = np.ndarray(shape=(N, N), dtype=float, order='F')
    X = np.ndarray(shape=(N), dtype=float, order='F')

    for i in range(N): 
        for j in range(N): A[i][j] = 0.0

    for i in range(N):
        A[i][i] = P2
        X[i] = 0
        if i != 0:
            A[i][i - 1] = P1
        if i != N-1:
            A[i][i + 1] = P1 

    X[N-1] = 1

    A, X = ml.gauss_elimination(A, X)
    X = X * (N * (N + 1) * (-1))
    print(X)

_N = int(sys.argv[1])
test_mm(_N)
test_gaussian(_N)