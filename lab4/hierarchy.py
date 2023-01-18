from typing import Union, List, Tuple
import numpy as np
from sklearn.utils.extmath import randomized_svd
from sys import getsizeof


# generate n x m matrix with 100*p % of values coming from uniform distribution on (0,1)
# remaining values should be initialized as zeros
def generate_matrix(size, fill):
    """Return matrix of size 2^n x 2^n with fill % of values coming from uniform distribution on (0,1)"""
    rand_matrix = np.random.rand(size, size)
    mask = np.random.choice([1, 0], size=(size, size), p=[fill, 1 - fill])
    return mask * rand_matrix


MatrixNodeData = Union[List[np.ndarray], List['MatrixNode']]


class MatrixNode:
    def __init__(self, is_leaf: bool, data: MatrixNodeData):
        # if true, the node contains compressed sub-matrix. otherwise reference to children nodes
        self.is_leaf = is_leaf

        # compressed sub-matrix or references
        self.data = data

    # d (size of decompressed d x d matrix)
    def get_size(self):
        if self.is_leaf:
            return len(self.data[0])
        else:
            return max(self.data[0][0].get_size(),
                       self.data[1][0].get_size(),
                       self.data[0][1].get_size(),
                       self.data[1][1].get_size()) * 2

    
    


def compress(matrix: np.ndarray, leaf_order=2, error=1e-5) -> MatrixNode:
    # d (size of the d x d input matrix)
    d = len(matrix)

    # recursive function used for compression takes coordinates of the sub-matrix in the input matrix
    def rek(dims: Tuple[int, int, int, int]) -> MatrixNode:

        # we extract the coordinates
        y0, y1, x0, x1 = dims

        # if it's 1 x 1 matrix then we should create a leaf node
        if y1 - y0 == 1:
            return MatrixNode(True, [np.array([matrix[y0, x0]]), np.array([[1]])])

        # we compress with truncated svd
        U, D, V_ = randomized_svd(matrix[y0:y1, x0:x1],
                                  n_components=leaf_order,
                                  random_state=42)

        # condition for compressing
        if D[-1] < error:
            # if satisfied we create a leaf node
            return MatrixNode(True, [U, np.diag(D) @ V_])

        else:
            # otherwise we need to evaluate the recursive function on sub-matrices
            data = [[],
                    []]

            data[0].append(rek((y0, (y1 + y0) // 2, x0, (x1 + x0) // 2)))
            data[0].append(rek((y0, (y1 + y0) // 2, (x1 + x0) // 2, x1)))
            data[1].append(rek(((y1 + y0) // 2, y1, x0, (x1 + x0) // 2)))
            data[1].append(rek(((y1 + y0) // 2, y1, (x1 + x0) // 2, x1)))

            # and create an internal node
            return MatrixNode(False, data)

    # we only need to call the recursive function on the full matrix
    return rek((0, d, 0, d))


def decompress(node: MatrixNode) -> np.ndarray:
    # we get the size of compressed matrix
    d = node.get_size()

    # initialize the result
    result = np.empty((d, d))

    # recursive function taking a node and coordinates of sub-matrix in the output matrix
    def rek(m_node, dims=(0, d, 0, d)):

        # we extract the coordinates
        y0, y1, x0, x1 = dims

        # if the node is a leaf we can multiply the elements of "data" field and copy to the result
        if m_node.is_leaf:
            result[y0:y1, x0:x1] = m_node.data[0] @ m_node.data[1]
        else:
            # otherwise we need to call the recursive function on the children nodes
            rek(m_node.data[0][0], (y0, (y1 + y0) // 2, x0, (x1 + x0) // 2))
            rek(m_node.data[0][1], (y0, (y1 + y0) // 2, (x1 + x0) // 2, x1))
            rek(m_node.data[1][0], ((y1 + y0) // 2, y1, x0, (x1 + x0) // 2))
            rek(m_node.data[1][1], ((y1 + y0) // 2, y1, (x1 + x0) // 2, x1))

    rek(node)
    return result


def demo() -> None:
    m = generate_matrix(2 ** 7, 0.2)
    mc = compress(m, 3)
    m1 = decompress(mc)
    print(((m - m1) ** 2).sum())


if __name__ == '__main__':
    demo()
