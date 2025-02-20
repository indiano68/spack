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
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=200)
# print(A)
# print()
# print(y)
# print()
# print(x)
# print(A@x)
# solution = np.linalg.solve(A,y)
# _, S, _ = np.linalg.svd(A)
# print(S[0]/S[-1])
# for idx, element in enumerate(x):
#     print(f"x[{idx}] = {element} {solution[idx]}")