import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from common.LinearCodes import LinearCodes

class PEG(LinearCodes):
    def __init__(self, rows, cols, l):
        self.n = cols
        self.m = rows
        self.l = l
        self.H = np.zeros((self.m, self.n), dtype=np.int8)
        self.D_s = np.zeros(self.n, dtype=np.int16)
        self.D_c = np.zeros(self.m, dtype=np.int16)
        self.D_s_future = np.zeros(self.n, dtype=np.int16) + l

    def _bfs(self, s):
        s_list = np.zeros(self.n, dtype = np.int8)
        s_list[s] = 1
        cur_c_list = np.zeros(self.m, dtype=np.int8)
        new_c_list = np.zeros(self.m, dtype=np.int8)
        c_Q = np.array((), dtype=np.int16)
        s_Q = np.array([s], dtype=np.int16)
        while(True):
            is_cycle = True
            if s_Q.size:
                for _s in s_Q:
                    cur_c_list_upd = np.argwhere(self.H[:, _s] == cur_c_list + 1).ravel()
                    new_c_list[self.H[:, _s] == 1] = 1
                    if cur_c_list_upd.size:
                        c_Q = np.hstack((c_Q, cur_c_list_upd)).astype(np.int16)
                        is_cycle = False
            s_Q = np.array(())
            if c_Q.size:
                for _c in c_Q:
                    s_list_update = np.argwhere(self.H[_c, :] == s_list + 1).ravel()
                    s_list[self.H[_c, :] == 1] = 1
                    if s_list_update.size:
                        s_Q = np.hstack((s_Q, s_list_update)).astype(np.int16)
            c_Q = np.array(())
            if is_cycle or np.count_nonzero(new_c_list) == self.m:
                ind = np.argwhere(cur_c_list == 0).ravel()
                degree = self.D_c[ind]
                return ind[np.argmin(degree)]
            cur_c_list = new_c_list.copy()

    def peg(self):
        for j in range(self.n):
            for k in range(self.D_s_future[j]):
                if k == 0:
                    c = np.argmin(self.D_c)
                else:
                    c = self._bfs(j)
                self.H[c, j] = 1
                self.D_s[j] += 1
                self.D_c[c] += 1
        np.savetxt(f'data/H_matrix/H{self.l}_{self.m}x{self.n}.txt', self.H, fmt='%d', delimiter=',')
