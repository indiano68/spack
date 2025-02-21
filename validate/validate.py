import numpy as np
import os
script_path = os.path.dirname(os.path.abspath(__file__))
prefix = script_path+'/gemv/test_gemv_'
suffix = '.txt'
A = np.loadtxt(prefix + 'A' + suffix)
x = np.loadtxt(prefix + 'x' + suffix)   
y = np.loadtxt(prefix + 'y' + suffix)

result = np.allclose(A @ x, y)
result_str = "\033[92mTrue\033[0m" if result else "\033[91mFalse\033[0m"

print("@ Testing gemv with given sizes: \n" + \
    "    A.shape = \033[94m" + str(A.shape) + "\033[0m\n" + \
    "    x.shape = \033[94m" + str(x.shape[0]) + "\033[0m\n" + \
    "    y.shape = \033[94m" + str(y.shape[0]) + "\033[0m\n" + \
    "-> RESULT: " + result_str)

prefix = script_path+'/backsob/test_'
suffix = '.txt'
A = np.loadtxt(prefix + 'upper_triangular' + suffix)
x = np.loadtxt(prefix + 'x' + suffix)   
y = np.loadtxt(prefix + 'b' + suffix)

result = np.allclose(np.linalg.solve(A,y), x)
result_str = "\033[92mTrue\033[0m" if result else "\033[91mFalse\033[0m"

print("@ Testing back_sobstitution with given sizes: \n" + \
    "    A.shape = \033[94m" + str(A.shape) + "\033[0m\n" + \
    "    x.shape = \033[94m" + str(x.shape[0]) + "\033[0m\n" + \
    "    y.shape = \033[94m" + str(y.shape[0]) + "\033[0m\n" + \
    "-> RESULT: " + result_str)

prefix = script_path+'/'
suffix = '.txt'
A = np.loadtxt(prefix + 'poisson_large' + suffix)
L = np.loadtxt(prefix + 'L' + suffix)   

L_numpy = np.linalg.cholesky(A)
result = np.allclose(L_numpy, L)
result_str = "\033[92mTrue\033[0m" if result else "\033[91mFalse\033[0m"
print("@ Testing cholesky decomposition with given sizes: \n" + \
    "    A.shape = \033[94m" + str(A.shape) + "\033[0m\n" + \
    "    L.shape = \033[94m" + str(L.shape) + "\033[0m\n" + \
    "-> RESULT: " + result_str)
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=200)
