import numpy as np 

create_row = lambda Z, A, En : (np.array([A, -A**(2/3), -(A-2*Z)**2/A, -Z**2/A**(1/3), (1-A%2)*(-1)**Z/(A**(3/4))]), En*A)
def create_matrix(data_matrix): 
    if data_matrix.shape[1] > 3 : raise Exception("Too many values to unpack")
    if data_matrix.shape[1] < 3 : raise Exception("Too less values to unpack")
    if data_matrix.shape[0] < 5 : raise Exception("Too less values to solve for coefficients, Pls give details of exactly 5 elemnts.")
    if data_matrix.shape[0] > 5 : raise Exception("Too many values to solve for coefficients, Pls give details of exactly 5 elements")
    A = np.zeros((5,5)) ; b = np.zeros((5))
    for i in range(5):
        eqn, const = create_row(data_matrix[i, 0], data_matrix[i, 1], data_matrix[i, 2])
        A[i, :] = eqn
        b[i] = const
    return np.linalg.solve(A,b)

data_matrix = np.array([[4, 10, 6.498], 
                        [6, 14, 7.520], 
                        [8, 18, 7.767],
                        [12, 26, 8.334],
                        [14, 30, 8.521]])
coeff = create_matrix(data_matrix)
print(coeff)
