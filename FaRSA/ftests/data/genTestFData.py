import numpy as np
from scipy.sparse import coo_matrix

# create regression data
np.random.seed(417)
m, n = 10, 6
A = np.random.randn(m, n)
b = np.random.rand(m, 1)
for i in range(m):
    A[i, i % n:i % n+5] = 0
Asparse = coo_matrix(A)
with open("lsMatrix.txt", 'w') as f:
    f.write(f"{m} {n} {Asparse.nnz}\n")
    for idx in range(len(Asparse.row)):
        i = Asparse.row[idx]
        j = Asparse.col[idx]
        f.write(f"{i} {j} {A[i,j]:+3.3e}\n")
with open("lsLabel.txt", 'w') as f:
    f.write(f"{m}\n")
    for idx in range(len(b)):
        f.write(f"{b[idx,0]:+3.3e}\n")
