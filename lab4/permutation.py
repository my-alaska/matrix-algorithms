import numpy as np
import matplotlib.pyplot as plt

from hierarchy import generate_matrix, compress
from total_size import total_size
import time


def minimum_degree_transformation(input_matrix: np.ndarray) -> np.ndarray:
    E = list(zip(*np.nonzero(input_matrix)))
    graph = [{i} for i in range(input_matrix.shape[0])]
    for e in E:
        graph[e[0]].add(e[1])
        graph[e[1]].add(e[0])

    order = []
    for _ in range(len(graph)):
        v = -1
        min_deg = len(input_matrix)

        for j, s in enumerate(graph):
            if len(s) != 0 and len(s) < min_deg:
                v = j
                min_deg = len(s)

        order.append(v)
        for s in graph:
            if v in s:
                s.remove(v)
        graph[v].clear()

    inverted = [None for _ in range(len(order))]
    for i, o in enumerate(order):
        inverted[o] = i

    result = np.zeros(input_matrix.shape)
    for coord in E:
        result[inverted[coord[0]]][inverted[coord[1]]] = input_matrix[coord]

    return result


def demo():
    matrix = generate_matrix(2**3, 0.5)
    print(matrix)
    print()
    print(minimum_degree_transformation(matrix))


def report():
    compression_results = []
    for i in range(4, 10):
        for _ in range(10):
            start = time.time()
            matrix = generate_matrix(2**i, 0.1)
            permuted = minimum_degree_transformation(matrix)
            compressed1 = compress(matrix)
            compressed2 = compress(permuted)
            size1 = total_size(compressed1)
            size2 = total_size(compressed2)
            size3 = total_size(matrix)
            size4 = total_size(permuted)

            stop = time.time()
            result = f"{i} {size1} {size2} {size3} {size4} {stop - start}"
            print(result)
            compression_results.append(result)


# function to estimate best parameters a and b for x^a * b given vectors x and y
def estimate_parameters(x, y):
    x = np.log(x)
    y = np.log(y)
    A = np.vstack([x, np.ones(len(x))]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return a, np.exp(b)


def generate_plots():
    with open("lab4/result.txt", "r") as f:
        results = [line.split() for line in f.readlines()]
        results = [(int(i), int(s1), int(s2), float(t)) for i, s1, s2, t in results]

    x = [2**i for i, s1, s1, t in results]
    y = [s1 for i, s1, s1, t in results]
    z = [s2 for i, s1, s2, t in results]

    a, b = estimate_parameters(x, y)
    c, d = estimate_parameters(x, z)
    X = np.linspace(0, 2**9, 100)

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.scatter(x, y, label="skompresowana oryginalna macierz", marker="x")
    ax.plot(X, X**a * b, label=f"y = x^{a:.2f} * {b:.2f}")
    ax.legend()
    ax.set_title("Wykres zajętej pamięci w zależności od rozmiaru d macierzy d*d\n dla skompresowanej oryginalnej macierzy")
    ax.set_xlabel("d")
    ax.set_ylabel("zajęta pamięć [B]")
    fig.savefig("lab4/compression1.png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.scatter(x, z, label="skompresowana permutowana macierz", marker="+")
    ax.plot(X, X**c * d, label=f"y = x^{c:.2f} * {d:.2f}")
    ax.legend()
    ax.set_title("Wykres zajętej pamięci w zależności od rozmiaru d macierzy d*d dla skompresowanej permutowanej macierzy")
    ax.set_xlabel("d")
    ax.set_ylabel("zajęta pamięć [B]")
    fig.savefig("lab4/compression2.png")

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.scatter(x, y, label="oryginalna macierz", marker="x")
    ax.scatter(x, z, label="permutowana macierz", marker="+")
    ax.legend()
    ax.set_title("Porównanie rozmiaru skompresowanych macierzy oryginalnej i permutowanej")
    ax.set_xlabel("d")
    ax.set_ylabel("zajęta pamięć [B]")
    fig.savefig("lab4/compression3.png")





if __name__ == "__main__":
    # np.set_printoptions(linewidth=200, precision=4, floatmode='fixed')

    # report()
    generate_plots()
