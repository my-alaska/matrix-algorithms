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


# print(generate_matrix(10,10,0.5))


class MNode:
    def __init__(self, is_leaf, data):
        # if true, the node contains compressed submatrix. otherwise reference to children nodes
        self.is_leaf = is_leaf

        # compressed submatrix or references
        self.data = data

    # d (size of decompressed d x d matrix)
    def get_size(self):
        if self.is_leaf:
            return len(self.data[0])
        else:
            return max(self.data[0][0].get_size(),
                       self.data[1][0].get_size(),
                       self.data[0][1].get_size(),
                       self.data[1][1].get_size())*2


def compress(M, r=2, e=1e-5):

    # d (size of the d x d input matrix)
    d = len(M)

    # recursive function used for compression takes coordinates of the submatrix in the input matrix
    def rek(dims):

        # we extract the coordinates
        y0, y1, x0, x1 = dims

        # if it's 1 x 1 matrix then we should create a leaf node
        if y1-y0 == 1: return MNode(True,[np.array([M[y0,x0]]),np.array([[1]])])

        # we compress with truncated svd
        U, D, V_ = randomized_svd(M[y0:y1, x0:x1],
                                  n_components=r,
                                  random_state=42)

        # condition for compressing
        if D[-1] < e:
            # if satisfied we create a leaf node
            return MNode(True, [U,np.diag(D)@V_])

        else:
            # otherwise we need to evaluate the recursive function on submatrices
            data = [[None, None],
                    [None, None]]

            data[0][0] = rek((y0, (y1 + y0) // 2, x0, (x1 + x0) // 2))
            data[0][1] = rek((y0, (y1 + y0) // 2, (x1 + x0) // 2, x1))
            data[1][0] = rek(((y1 + y0) // 2, y1, x0, (x1 + x0) // 2))
            data[1][1] = rek(((y1 + y0) // 2, y1, (x1 + x0) // 2, x1))

            # and create an internal node
            return MNode(False, data)

    # we only need to call the recursive function on the full matrix
    return rek((0, d, 0, d))


def decompress(mnode):
    # we get the size of compressed matrix
    d = mnode.get_size()

    # initialize the result
    result = np.empty((d,d))

    # recursive function taking a node and coordinates of submatrix in the output matrix
    def rek(mnode,dims=(0,d,0,d)):

        # we extract the coordinates
        y0, y1, x0, x1 = dims

        # if the node is a leaf we can multiply the elements of "data" field and copy to the result
        if mnode.is_leaf:
            result[y0:y1, x0:x1] = mnode.data[0]@mnode.data[1]
        else:
            # otherwise we need to call the recursive function on the children nodes
            rek(mnode.data[0][0], (y0, (y1 + y0) // 2, x0, (x1 + x0) // 2))
            rek(mnode.data[0][1], (y0, (y1 + y0) // 2, (x1 + x0) // 2, x1))
            rek(mnode.data[1][0], ((y1 + y0) // 2, y1, x0, (x1 + x0) // 2))
            rek(mnode.data[1][1], ((y1 + y0) // 2, y1, (x1 + x0) // 2, x1))

    rek(mnode)
    return result

m = generate_matrix(2**7, 2**7, 0.2)
mc = compress(m,3)
m1 = decompress(mc)

print(((m-m1)**2).sum())