from math import log, ceil
import random
import time


def binet(A, B, d=None):
    if d == None:
        K, L, M = len(A), len(B), len(B[0])
    else:
        K, L, M = d

    if K == 0:
        return [], 0
    if L == 0:
        return [[0 for j in range(M)] for i in range(K)], 0
    if M == 0:
        return [[]], 0
    if K == 1 and L == 1 and M == 1:
        return [[A[0][0] * B[0][0]]], 1

    A11 = [[A[i][j] for j in range(0, L // 2)] for i in range(0, K // 2)]
    A12 = [[A[i][j] for j in range(L // 2, L)] for i in range(0, K // 2)]
    A21 = [[A[i][j] for j in range(0, L // 2)] for i in range(K // 2, K)]
    A22 = [[A[i][j] for j in range(L // 2, L)] for i in range(K // 2, K)]

    B11 = [[B[i][j] for j in range(0, M // 2)] for i in range(0, L // 2)]
    B12 = [[B[i][j] for j in range(M // 2, M)] for i in range(0, L // 2)]
    B21 = [[B[i][j] for j in range(0, M // 2)] for i in range(L // 2, L)]
    B22 = [[B[i][j] for j in range(M // 2, M)] for i in range(L // 2, L)]

    C = [[0 for j in range(M)] for i in range(K)]

    A11B11, c1 = binet(A11, B11, (K // 2, L // 2, M // 2))
    A12B21, c2 = binet(A12, B21, (K // 2, L - L // 2, M // 2))

    A11B12, c3 = binet(A11, B12, (K // 2, L // 2, M - M // 2))
    A12B22, c4 = binet(A12, B22, (K // 2, L - L // 2, M - M // 2))

    A21B11, c5 = binet(A21, B11, (K - K // 2, L // 2, M // 2))
    A22B21, c6 = binet(A22, B21, (K - K // 2, L - L // 2, M // 2))

    A21B12, c7 = binet(A21, B12, (K - K // 2, L // 2, M - M // 2))
    A22B22, c8 = binet(A22, B22, (K - K // 2, L - L // 2, M - M // 2))

    for i in range(K // 2):
        for j in range(M // 2):
            C[i][j] += A11B11[i][j] + A12B21[i][j]

    for i in range(K // 2):
        for j in range(M - M // 2):
            C[i][M // 2 + j] += A11B12[i][j] + A12B22[i][j]

    for i in range(K - K // 2):
        for j in range(M // 2):
            C[K // 2 + i][j] += A21B11[i][j] + A22B21[i][j]

    for i in range(K - K // 2):
        for j in range(M - M // 2):
            C[K // 2 + i][M // 2 + j] += A21B12[i][j] + A22B22[i][j]

    cost = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + K * M
    return C, cost






def add_m(A, B):
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError
    X, Y = len(A), len(A[0])
    return [[A[row][col]+B[row][col] for col in range(Y)] for row in range(X)]


def sub_m(A, B):
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError
    X, Y = len(A), len(A[0])
    return [[A[row][col]-B[row][col] for col in range(Y)] for row in range(X)]


def strassen_r(A, B):
    K, L, M = len(A), len(B), len(B[0])

    if K == 1 and L == 1 and M == 1:
        return [[A[0][0] * B[0][0]]], 1

    A11 = [[A[i][j] for j in range(0, L // 2)] for i in range(0, K // 2)]
    A12 = [[A[i][j] for j in range(L // 2, L)] for i in range(0, K // 2)]
    A21 = [[A[i][j] for j in range(0, L // 2)] for i in range(K // 2, K)]
    A22 = [[A[i][j] for j in range(L // 2, L)] for i in range(K // 2, K)]

    B11 = [[B[i][j] for j in range(0, M // 2)] for i in range(0, L // 2)]
    B12 = [[B[i][j] for j in range(M // 2, M)] for i in range(0, L // 2)]
    B21 = [[B[i][j] for j in range(0, M // 2)] for i in range(L // 2, L)]
    B22 = [[B[i][j] for j in range(M // 2, M)] for i in range(L // 2, L)]

    C = [[0 for j in range(M)] for i in range(K)]
    S = 0

    P1, c1 = strassen_r(add_m(A11, A22), add_m(B11, B22))
    S += c1 + (len(A11)*len(A11[0]) + len(B11)*len(B11[0])) # number of floating point additions for the poperation

    P2, c2 = strassen_r(add_m(A21, A22), B11)
    S += c2 + len(A21)*len(A21[0])

    P3, c3 = strassen_r(A11, sub_m(B12, B22))
    S += c3 + len(B12)*len(B12[0])

    P4, c4 = strassen_r(A22, sub_m(B21, B11))
    S += c4 + len(B12)*len(B12[0])

    P5, c5 = strassen_r(add_m(A11, A12), B22)
    S += len(A11)*len(A11[0])

    P6, c6 = strassen_r(sub_m(A21, A11), add_m(B11, B12))
    S += c6 + (len(A21) * len(A21[0]) + len(B11) * len(B11[0]))

    P7, c7 = strassen_r(sub_m(A12, A22), add_m(B21, B22))
    S += c7 + (len(A12)*len(A12[0]) + len(B21)*len(B21[0]))



    C11 = add_m(sub_m(add_m(P1, P4), P5), P7)
    S += 3 * len(P1)*len(P1[0])

    C12 = add_m(P3, P5)
    S += 3 * len(P3)*len(P3[0])

    C21 = add_m(P2, P4)
    S += 3 * len(P2) * len(P2[0])

    C22 = add_m(add_m(sub_m(P1, P2), P3), P6)
    S += 3 * len(P1) * len(P1[0])

    for i in range(K // 2):
        for j in range(M // 2):
            C[i][j] += C11[i][j]

    for i in range(K // 2):
        for j in range(M - M // 2):
            C[i][M//2+j] += C12[i][j]

    for i in range(K - K // 2):
        for j in range(M // 2):
            C[K//2+i][j] = C21[i][j]

    for i in range(K - K // 2):
        for j in range(M - M // 2):
            # C[K // 2 + i][M // 2 + j] += A21B12[i][j] + A22B22[i][j]
            C[K//2 + i][M//2+j] = C22[i][j]

    return C, S


def strassen(A, B):
    k = max(len(A), len(B), len(B[0]))
    k = 2**int(ceil(log(k, 2)))

    A0 = [[0 for _ in range(k)] for _ in range(k)]
    B0 = [[0 for _ in range(k)] for _ in range(k)]

    for i in range(len(A)):
        for j in range(len(A[0])):
            A0[i][j] = A[i][j]

    for i in range(len(B)):
        for j in range(len(B[0])):
            B0[i][j] = B[i][j]

    C0, S = strassen_r(A0, B0)
    C = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            C[i][j] = C0[i][j]

    return C, S


def print_arr(A):
    for row in A:
        print(row)


def test():
    A = [[1, 1, 1, 1],
         [2, 2, 2, 2],
         [3, 3, 3, 3]]

    B = [[1, 1],
         [2, 2],
         [3, 3],
         [4, 4]]

    C = [[1, 2],
         [3, 4]]

    D = [[1, 2],
         [1, 2]]

    E = [[1, 2, 3, 4],
         [1, 2, 3, 4],
         [1, 2, 3, 4],
         [1, 2, 3, 4]]

    F = [[1, 1, 1, 1],
         [2, 2, 2, 2],
         [3, 3, 3, 3],
         [4, 4, 4, 4]]

    a, ac = binet(C, D)
    b, bc = strassen_r(C, D)
    c, cc = strassen(C, D)
    assert a == b == c

    a, ac = binet(E, F)
    b, bc = strassen_r(E, F)
    c, cc = strassen(E, F)
    assert a == b == c

    a, ac = binet(A, B)
    c, cc = strassen(A, B)
    assert a == c


def measure():
    # file_strassen = open('strassen.txt', 'w')
    # file_binet = open('binet.txt', 'w')

    random.seed(42)
    for i in range(2, 10):
        k = 2**i
        A = [[random.random() for n in range(k)] for m in range(k)]
        B = [[random.random() for n in range(k)] for m in range(k)]

        start = time.time()
        C, dc = strassen(A, B)
        stop = time.time()
        # file_strassen.write(f'i: {i} {stop-start}\n')
        print(f'strassen, i: {i}, k: {k}, time: {stop-start}, operations: {dc}')

        start = time.time()
        D, dc = binet(A, B)
        stop = time.time()

        # file_binet.write(f'i: {i} {stop-start}\n')
        print(f'binet, i: {i}, k: {k}, time: {stop-start}, operations: {dc}')

    # file_strassen.close()
    # file_binet.close()


if __name__ == '__main__':
    test()
    measure()
