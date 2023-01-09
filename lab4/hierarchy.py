import numpy as np
from sklearn.utils.extmath import randomized_svd


# generate n x m matrix with 100*p % of values coming from uniform distribution on (0,1)
# remaining values should be initialized as zeros
def generate_matrix(y, x, p):
    # random values between 0 and 1
    rand_matrix = np.random.rand(y, x)

    # mask consisting of 100*p % of zeros at random positions
    mask = np.random.choice([1, 0], size=(y, x), p=[p, 1 - p])

    # element-wise multiplication
    return mask * rand_matrix


def minimum_degree(M):

    E = list(zip(*np.nonzero(M)))
    graph = [{i} for i in range(M.shape[0])]
    for e in E:
        graph[e[0]].add(e[1])
        graph[e[1]].add(e[0])


    order = []
    for _ in range(len(graph)):
        v = -1
        min_deg = len(M)

        for j,s in enumerate(graph):
            if len(s) != 0 and len(s) < min_deg:
                v = j
                min_deg = len(s)

        order.append(v)
        for s in graph:
            if v in s:
                s.remove(v)
        graph[v].clear()

    # print(order)

    inverted = [None for _ in range(len(order))]
    for i,o in enumerate(order):
        inverted[o] = i

    result = np.zeros(M.shape)
    for coord in E:
        result[inverted[coord[0]]][inverted[coord[1]]] = M[coord]

    return result




m = generate_matrix(2 ** 3, 2 ** 3, 0.1)
print(m)
print()
print(minimum_degree(m))