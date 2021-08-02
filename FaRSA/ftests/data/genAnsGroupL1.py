from scipy.sparse import csr_matrix
import numpy as np


class GL1:
    def __init__(self, Lambda, group):
        """
        !!Warning: need `group` be ordered in a consecutive manner, i.e., 
        group: array([1., 1., 1., 2., 2., 2., 3., 3., 3., 3.])
        Then:
        unique_groups: array([1., 2., 3.])
        group_frequency: array([3, 3, 4]))
        """
        self.group = group
        self.Lambda = Lambda
        self.unique_groups, self.group_frequency = np.unique(
            self.group, return_counts=True)
        self.Lambda_group = self.Lambda * np.sqrt(self.group_frequency)
        self.K = len(self.unique_groups)
        self.group_size = -1 * np.ones(self.K)
        p = group.shape[0]
        full_index = np.arange(p)
        starts = []
        ends = []
        data = np.zeros(p)
        row_idx = np.zeros(p)
        col_idx = np.zeros(p)
        for i in range(self.K):
            G_i = full_index[np.where(self.group == self.unique_groups[i])]
            # record the `start` and `end` indices of the group G_i to avoid fancy indexing innumpy
            # in the example above, the start index and end index for G_1 is 0 and 2 respectively
            # since python `start:end` will include `start` and exclude `end`, so we will add 1 to the `end`
            # so the G_i-th block of X is indexed by X[start:end]
            start, end = min(G_i), max(G_i) + 1
            starts.append(start)
            ends.append(end)
            self.group_size[i] = end - start
            data[start:end] = self.Lambda_group[i] ** 2
            row_idx[start:end] = i
            col_idx[start:end] = np.arange(start, end)
        # wrap as np.array for jit compile purpose
        self.starts = np.array(starts)
        self.ends = np.array(ends)
        self.mat = csr_matrix((data, (row_idx, col_idx)), shape=(self.K, p))

    def __str__(self):
        return("Group L1")

    def func(self, X):
        return np.sum(np.sqrt(self.mat@(X ** 2)))

    def grad(self, X, sub_grp_idx):
        g = []
        for i in sub_grp_idx:
            start, end = self.starts[i], self.ends[i]
            X_Gi = X[start:end]
            X_Gi_norm = np.sqrt(np.sum(X_Gi*X_Gi))
            g_Gi = self.Lambda_group[i] * X_Gi / X_Gi_norm
            for e in g_Gi:
                g.append(e)
        return np.array(g)

    def proximal(self, x, gradfx, stepsize):
        prox = np.zeros_like(x)
        gradstep = x - stepsize * gradfx
        for i in range(self.K):
            start, end = self.starts[i], self.ends[i]
            gradstep_Gi = gradstep[start:end]
            prox[start:end] = max(0, 1 - (self.Lambda_group[i] * stepsize/np.sqrt(np.sum(gradstep_Gi*gradstep_Gi))
                                          )
                                  ) * gradstep_Gi
        return prox

    def _prepare_hv_data(self, X, subgroup_index):
        self.hv_data = {}
        start = 0
        for i in subgroup_index:
            start_x, end_x = self.starts[i], self.ends[i]
            XG_i = X[start_x:end_x]
            XG_i_norm = np.sqrt(np.dot(XG_i.T, XG_i))[0][0]
            end = start + end_x - start_x
            self.hv_data[i] = {}
            self.hv_data[i]['XG_i'] = XG_i
            self.hv_data[i]['XG_i_norm'] = XG_i_norm
            self.hv_data[i]['start'] = start
            self.hv_data[i]['end'] = end
            self.hv_data[i]['XG_i_norm_cubic'] = XG_i_norm**3
            start = end

    def hessian_vector_product_fast(self, v, subgroup_index):
        hv = np.empty_like(v)
        for i in subgroup_index:
            start = self.hv_data[i]['start']
            end = self.hv_data[i]['end']
            vi = v[start:end]
            temp = np.matmul(self.hv_data[i]['XG_i'].T, vi)
            hv[start:end] = self.Lambda_group[i] * (1 / self.hv_data[i]['XG_i_norm'] * vi -
                                                    (temp / self.hv_data[i]['XG_i_norm_cubic']) *
                                                    self.hv_data[i]['XG_i'])
        return hv


r = GL1(1.23, np.array([0, 0, 1, 1, 1, 2]))
x = np.array([0.0, 0.0, 1.1, 0.0, 2.2, 3.3]).reshape(-1, 1)
rfun = r.func(x)
rgrad = r.grad(x, [1, 2])
rprox1 = r.proximal(x, x + 0.1, 0.2)
rprox2 = r.proximal(x, x + 2.0, 0.2)
cols = [2, 3, 4, 5]
v = np.array([i*1.1 for i in range(len(cols))]).reshape(-1, 1)
r._prepare_hv_data(x, [1, 2])
rHv = r.hessian_vector_product_fast(v, [1, 2])

with open("rfun.txt", 'w') as f:
    f.write(f"{1}\n")
    f.write(f"{rfun:+3.9e}\n")
with open("rgrad.txt", 'w') as f:
    f.write(f"{rgrad.shape[0]}\n")
    for idx in range(len(rgrad)):
        f.write(f"{rgrad[idx,0]:+3.9e}\n")
with open("rprox1.txt", 'w') as f:
    f.write(f"{rprox1.shape[0]}\n")
    for idx in range(len(rprox1)):
        f.write(f"{rprox1[idx,0]:+3.9e}\n")
with open("rprox2.txt", 'w') as f:
    f.write(f"{rprox2.shape[0]}\n")
    for idx in range(len(rprox2)):
        f.write(f"{rprox2[idx,0]:+3.9e}\n")
with open("rHv.txt", 'w') as f:
    f.write(f"{rHv.shape[0]}\n")
    for idx in range(len(rHv)):
        f.write(f"{rHv[idx,0]:+3.9e}\n")
