from lab1.matrix_multiplication import strassen
import numpy as np
import time

def inv(M):
    k = len(M)//2
    if k == 0:
        return [[1/M[0][0]]], 1

    # splitting the input into A11, A12, A21, A22
    A11 = [[ M[i][j] for j in range(k)] for i in range(k)]
    A12 = [[ M[i][j] for j in range(k,2*k)] for i in range(k)]
    A21 = [[ M[i][j] for j in range(k)] for i in range(k,2*k)]
    A22 = [[ M[i][j] for j in range(k,2*k)] for i in range(k,2*k)]

    # inverting A11 and initializing p
    invA11, p = inv(A11)

    # calculating S22
    Q1, i = strassen(invA11,A12)
    p += i
    Q2, i = strassen(A21,Q1)
    p += i
    S22 = [[ A22[i][j] - Q2[i][j] for j in range(k)] for i in range(k)]
    p += k**2

    # inverting S22
    invS22, i = inv(S22)
    p += i

    # getting the top left result
    Q1, i = strassen(A21,invA11)
    p += i
    Q2, i = strassen(invS22,Q1)
    p += i
    Q3, i = strassen(A12,Q2)
    p += i
    Q4, i = strassen(invA11,Q3)
    p += i
    B11 = [[ invA11[i][j] + Q4[i][j] for j in range(k)] for i in range(k)]
    p += k**2

    # getting the negated top right result
    Q1, i = strassen(A12,invS22)
    p += i
    negB12, i = strassen(invA11,Q1)
    p += i

    # getting the negated bottom left result
    Q1, i = strassen(A21,invA11)
    p +=i
    negB21, i = strassen(invS22,Q1)
    p+= i


    # initializing the result
    B = [[None for i in range(2*k)]for j in range(2*k)]

    # bottom right corner
    for i in range(k):
        for j in range(k):
            B[i][j] = B11[i][j]
            B[k+i][k+j] = invS22[i][j]
            B[i][k+j] = -negB12[i][j]
            B[k+i][j] = -negB21[i][j]

    return B, p


def lu(M):
    k = len(M) // 2
    if k == 0:
        return [[1]], M, 0

    # splitting the input into A11, A12, A21, A22
    A11 = [[M[i][j] for j in range(k)] for i in range(k)]
    A12 = [[M[i][j] for j in range(k, 2 * k)] for i in range(k)]
    A21 = [[M[i][j] for j in range(k)] for i in range(k, 2 * k)]
    S = [[M[i][j] for j in range(k, 2 * k)] for i in range(k, 2 * k)]

    # factorizing A11 into L11 and U11
    L11, U11, p = lu(A11)

    # inverting L11 and U11
    invL11,i = inv(L11)
    p += i
    invU11,i = inv(U11)
    p += i

    # calculating U12 and L21
    U12, i = strassen(invL11,A12)
    p+=i
    L21, i = strassen(A21,invU11)
    p+=i

    # calculating S
    Q1, i = strassen(L21,U12)
    p+=i
    for i in range(k):
        for j in range(k):
            S[i][j] -= Q1[i][j]
    p+= k**2

    # factorizing S into Ls and Us
    Ls, Us, i = lu(S)
    p+=i

    # initializing and filling L and U (result)
    L = [[0 for i in range(2*k)] for i in range(2*k)]
    U = [[0 for i in range(2 * k)] for i in range(2 * k)]
    for i in range(k):
        for j in range(k):
            L[i][j] = L11[i][j]
            L[k + i][k + j] = Ls[i][j]
            L[k + i][j] = L21[i][j]
            U[i][j] = U11[i][j]
            U[k+i][k+j] = Us[i][j]
            U[i][k+j] = U12[i][j]
    return L, U, p


def det(M):
    L, U, p = lu(M)
    n = len(M)
    result = U[0][0]
    for i in range(1, n):
        result *= U[i][i]
    return result, p + n - 1


def test():
    matrix = [[1,-2,-2,-3],
            [3,-9,0,-9],
            [-1,2,4,7],
            [-3,-6,26,2]]
    L, U, p = lu(matrix)
    # for i in range(len(L)):
    #     print(L[i], "\t\t", U[i])

    assert np.isclose(np.array([[1, 0, 0, 0], [3, 1, 0, 0], [-1, 0, 1, 0], [-3, 4, -2, 1]]), np.array(L)).all()
    assert np.isclose(np.array([[1, -2, -2, -3], [0, -3, 6, 0], [0, 0, 2, 4], [0, 0, 0, 1]]), np.array(U)).all()

def measure_lu(loops=10):
    for i in range(0, loops):
        k = 2**i
        while True:
            A = [[np.random.rand() for n in range(k)] for m in range(k)]
            B = np.array(A)
            if np.linalg.matrix_rank(B) == len(B):
                break
        
        start = time.time()
        L, U, p = lu(A)
        stop = time.time()

        print(f'{k} {stop-start} {p}')


def measure_inv(loops=10):
        for i in range(0, loops):
            k = 2**i
            while True:
                A = [[np.random.rand() for n in range(k)] for m in range(k)]
                B = np.array(A)
                if np.linalg.matrix_rank(B) == len(B):
                    break
            
            start = time.time()
            C, p = inv(A)
            stop = time.time()

            print(f'{k} {stop-start} {p}')


def measure_det(loops=10):
        for i in range(0, loops):
            k = 2**i
            while True:
                A = [[np.random.rand() for n in range(k)] for m in range(k)]
                B = np.array(A)
                if np.linalg.matrix_rank(B) == len(B):
                    break
            
            start = time.time()
            # print(A)
            C, p = det(A)
            
            stop = time.time()

            print(f'{k} {stop-start} {p}')


if __name__ == '__main__':
    test()
    measure_lu(8)
