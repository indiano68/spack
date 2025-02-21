import numpy as np

def build_five_point_stencil(n):
    N = n * n
    A = np.zeros((N, N), dtype=float)
    indices = np.arange(N).reshape(n, n)
    
    A[indices, indices] = 4.0
    
    # Left and right neighbors
    A[indices[:, :-1], indices[:, 1:]] = -1.0
    A[indices[:, 1:], indices[:, :-1]] = -1.0
    
    # Top and bottom neighbors
    A[indices[:-1, :], indices[1:, :]] = -1.0
    A[indices[1:, :], indices[:-1, :]] = -1.0
    
    # Apply Dirichlet boundary conditions
    # for i in range(n):
    #     A[indices[i, 0], :] = 0
    #     A[indices[i, 0], indices[i, 0]] = 1
    #     A[indices[i, -1], :] = 0
    #     A[indices[i, -1], indices[i, -1]] = 1
    #     A[indices[0, i], :] = 0
    #     A[indices[0, i], indices[0, i]] = 1
    #     A[indices[-1, i], :] = 0
    #     A[indices[-1, i], indices[-1, i]] = 1
    
    return A

np.set_printoptions(threshold=np.inf, linewidth=200)
matrix = build_five_point_stencil(4)
print(matrix)
np.set_printoptions(threshold=np.inf,precision =2, linewidth=200)
chle = np.linalg.cholesky(matrix)
print(chle)

# Write the matrix to a file
with open('matrix_output.txt', 'w') as f:
    for row in matrix:
        np.savetxt(f, row[np.newaxis], fmt='%.2f')