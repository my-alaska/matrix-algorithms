from lab1.matrix_multiplication import strassen

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

matrix = [[1,0,0,0],
          [0,3,0,0],
          [0,0,1,0],
          [0,0,0,2]]

print(inv(matrix))
