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

# TODO Strassen's methode

# A = [[1, 1, 1, 1],
#      [2, 2, 2, 2],
#      [3, 3, 3, 3]]
#
# B = [[1, 1],
#      [2, 2],
#      [3, 3],
#      [4, 4]]
#
# print(binet(A, B))


def strassen(A, B):
    K, L, M = len(A), len(B), len(B[0])
    # TODO


def print_arr(A):
    for row in A:
        print(row)


if __name__ == '__main__':

    A = [[1, 1, 1, 1],
         [2, 2, 2, 2],
         [3, 3, 3, 3]]

    B = [[1, 1],
         [2, 2],
         [3, 3],
         [4, 4]]

    a = binet(A, B)
    print_arr(a)
