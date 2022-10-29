from math import log, ceil


def binet(A, B):
    K, L, M = len(A), len(B), len(B[0])

    if K == 0 or L == 0 or M == 0:
        return [[]]
    if K == 1 and L == 1 and M == 1:
        return [[A[0][0] * B[0][0]]]

    A11 = [[A[i][j] for j in range(0, L // 2)] for i in range(0, K // 2)]
    A12 = [[A[i][j] for j in range(L // 2, L)] for i in range(0, K // 2)]
    A21 = [[A[i][j] for j in range(0, L // 2)] for i in range(K // 2, K)]
    A22 = [[A[i][j] for j in range(L // 2, L)] for i in range(K // 2, K)]

    B11 = [[B[i][j] for j in range(0, M // 2)] for i in range(0, L // 2)]
    B12 = [[B[i][j] for j in range(M // 2, M)] for i in range(0, L // 2)]
    B21 = [[B[i][j] for j in range(0, M // 2)] for i in range(L // 2, L)]
    B22 = [[B[i][j] for j in range(M // 2, M)] for i in range(L // 2, L)]

    C = [[0 for j in range(M)] for i in range(K)]

    A11B11 = binet(A11, B11)
    A12B21 = binet(A12, B21)

    A11B12 = binet(A11, B12)
    A12B22 = binet(A12, B22)

    A21B11 = binet(A21, B11)
    A22B21 = binet(A22, B21)

    A21B12 = binet(A21, B12)
    A22B22 = binet(A22, B22)

    for i in range(K // 2):
        for j in range(M // 2):
            C[i][j] += A11B11[i][j] + A12B21[i][j]

    for i in range(K // 2):
        for j in range(M - M // 2):
            C[i][M//2+j] += A11B12[i][j] + A12B22[i][j]

    for i in range(K - K // 2):
        for j in range(M // 2):
            C[K // 2 + i][j] += A21B11[i][j] + A22B21[i][j]

    for i in range(K - K // 2):
        for j in range(M - M // 2):
            C[K // 2 + i][M // 2 + j] += A21B12[i][j] + A22B22[i][j]

    return C


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
        return [[A[0][0] * B[0][0]]]

    A11 = [[A[i][j] for j in range(0, L // 2)] for i in range(0, K // 2)]
    A12 = [[A[i][j] for j in range(L // 2, L)] for i in range(0, K // 2)]
    A21 = [[A[i][j] for j in range(0, L // 2)] for i in range(K // 2, K)]
    A22 = [[A[i][j] for j in range(L // 2, L)] for i in range(K // 2, K)]

    B11 = [[B[i][j] for j in range(0, M // 2)] for i in range(0, L // 2)]
    B12 = [[B[i][j] for j in range(M // 2, M)] for i in range(0, L // 2)]
    B21 = [[B[i][j] for j in range(0, M // 2)] for i in range(L // 2, L)]
    B22 = [[B[i][j] for j in range(M // 2, M)] for i in range(L // 2, L)]

    C = [[0 for j in range(M)] for i in range(K)]

    P1 = strassen_r(add_m(A11, A22), add_m(B11, B22))
    P2 = strassen_r(add_m(A21, A22), B11)
    P3 = strassen_r(A11, sub_m(B12, B22))
    P4 = strassen_r(A22, sub_m(B21, B11))
    P5 = strassen_r(add_m(A11, A12), B22)
    P6 = strassen_r(sub_m(A21, A11), add_m(B11, B12))
    P7 = strassen_r(sub_m(A12, A22), add_m(B21, B22))

    C11 = add_m(sub_m(add_m(P1, P4), P5), P7)
    C12 = add_m(P3, P5)
    C21 = add_m(P2, P4)
    C22 = add_m(add_m(sub_m(P1, P2), P3), P6)

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

    return C


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

    C0 = strassen_r(A0, B0)
    C = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            C[i][j] = C0[i][j]

    return C


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

    a = binet(C, D)
    b = strassen_r(C, D)
    c = strassen(C, D)
    assert a == b == c

    a = binet(E, F)
    b = strassen_r(E, F)
    c = strassen(E, F)
    assert a == b == c

    a = binet(A, B)
    c = strassen(A, B)
    assert a == c


if __name__ == '__main__':
    test()
