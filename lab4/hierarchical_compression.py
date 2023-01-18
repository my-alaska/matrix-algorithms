import numpy as np
from sklearn.utils.extmath import randomized_svd
from sys import getsizeof
from typing import Tuple
import matplotlib.pyplot as plt


def generate_matrix(size, fill):
    """Return matrix of size 2^n x 2^n with fill % of values coming from uniform distribution on (0,1)"""
    rand_matrix = np.random.rand(size, size)
    mask = np.random.choice([1, 0], size=(size, size), p=[fill, 1 - fill])
    return mask * rand_matrix


class Node:
    pass


class NonCompressedLeafNode(Node):
    def __init__(self, matrix):
        self.matrix = matrix
        self.shape = matrix.shape

    def print(self, indent):
        print(f"{'  ' * indent}Leaf node of size {self.matrix.shape}")
        for row in self.matrix:
            print(f"{'  ' * (indent + 1)}{row}")

    def __sizeof__(self):
        return self.matrix.nbytes


class CompressedLeafNode(Node):
    def __init__(self, u: np.ndarray, v: np.ndarray):
        self.u = u
        self.v = v
        self.shape = u.shape[0], v.shape[1]

    def print(self, indent):
        print(f"{'  ' * indent}Compressed leaf node")
        print(f"{'  ' * (indent + 1)}U:")
        for row in self.u:
            print(f"{'  ' * (indent + 1)}{row}")
        print(f"{'  ' * (indent + 1)}V:")
        for row in self.v:
            print(f"{'  ' * (indent + 1)}{row}")

    def __sizeof__(self):
        return self.u.nbytes + self.v.nbytes


class CompressedInternalNode(Node):
    def __init__(self, a: Node, b: Node, c: Node, d: Node, shape: Tuple[int, int]):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.shape = shape

    def print(self, indent):
        print(f"{'  ' * indent}Internal node")
        self.a.print(indent + 1)
        self.b.print(indent + 1)
        self.c.print(indent + 1)
        self.d.print(indent + 1)

    def __sizeof__(self):
        return (
            getsizeof(self.a)
            + getsizeof(self.b)
            + getsizeof(self.c)
            + getsizeof(self.d)
        )


A = np.array(
    [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]], dtype=np.float64
)
B = np.array([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], dtype=np.float64)


def compress(
    matrix: np.ndarray, k: int = 2, epsilon: float = 1e-12, depth=0
) -> CompressedLeafNode:
    d = len(matrix)
    if k > d - 1:
        return NonCompressedLeafNode(matrix.copy())

    u, s, v = randomized_svd(matrix, n_components=k + 1, n_iter=5, random_state=42)
    u = u[:, :-1]
    s = s[:-1]
    v = v[:-1, :]
    if s[-1] < epsilon:
        return CompressedLeafNode(u, np.diag(s) @ v)
    else:
        a = compress(matrix[0 : d // 2, 0 : d // 2], k, epsilon, depth=depth + 1)
        b = compress(matrix[0 : d // 2, d // 2 : d], k, epsilon, depth=depth + 1)
        c = compress(matrix[d // 2 : d, 0 : d // 2], k, epsilon, depth=depth + 1)
        d = compress(matrix[d // 2 : d, d // 2 : d], k, epsilon, depth=depth + 1)
        return CompressedInternalNode(a, b, c, d, matrix.shape)


def decompress(node: Node, matrix: np.ndarray, coords: Tuple[int, int]) -> np.ndarray:
    if isinstance(node, NonCompressedLeafNode):
        matrix[
            coords[0] : coords[0] + node.matrix.shape[0],
            coords[1] : coords[1] + node.matrix.shape[1],
        ] = node.matrix
    elif isinstance(node, CompressedLeafNode):
        matrix[
            coords[0] : coords[0] + node.u.shape[0],
            coords[1] : coords[1] + node.v.shape[1],
        ] = (
            node.u @ node.v
        )
    elif isinstance(node, CompressedInternalNode):
        decompress(node.a, matrix, (coords[0], coords[1]))
        decompress(node.b, matrix, (coords[0], coords[1] + node.b.shape[1] // 2))
        decompress(node.c, matrix, (coords[0] + node.c.shape[0] // 2, coords[1]))
        decompress(
            node.d,
            matrix,
            (coords[0] + node.d.shape[0] // 2, coords[1] + node.d.shape[1] // 2),
        )

# write function node_representation that takes as input Node and returns a matrix that represents the node
# if the node is a compressed leaf node, write the matrix as matrix full of zeros, but first k rows and columns are filled with ones.
# k is the length of vector v in the compressed leaf node 
# do the same if the node is a compressed zero node. If node is a non compressed leaf node, write the matrix as it is.
# if node is an internal node, write the matrix as stack of matrices of its children

def sparsity_representation(node: Node) -> np.ndarray:
    if isinstance(node, NonCompressedLeafNode):
        return np.ones(node.matrix.shape)
    elif isinstance(node, CompressedLeafNode):
        k = node.v.shape[0]
        matrix = np.zeros(node.shape)
        matrix[:k, :] = 1
        matrix[:, :k] = 1
        return matrix
    elif isinstance(node, CompressedInternalNode):
        a = sparsity_representation(node.a)
        b = sparsity_representation(node.b)
        c = sparsity_representation(node.c)
        d = sparsity_representation(node.d)
        return np.vstack((np.hstack((a, b)), np.hstack((c, d))))



print("Here")
print()
A = generate_matrix(32, 0.1)
print(A)
print("Here 2")
B = compress(A, k=5)
B.print(0)

print("Here 3")
C = decompress(B, np.zeros(A.shape), (0, 0))

print(sparsity_representation(B))
plt.imshow(sparsity_representation(B), cmap="gray_r")
plt.show()