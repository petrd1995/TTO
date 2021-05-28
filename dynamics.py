from scipy import linalg as la
import numpy as np


def eigen(A, B):
    '''
        Computation and sorting of eigenvalues and eigenvectors from the smallest to the biggest
    '''

    eigenValues, eigenVectors = la.eig(A, B)
    idx = np.argsort(eigenValues)
    # odmocnina kvůli tomu, že omega (vl.frekvence) je umocněná
    eigenValues = np.around(np.sqrt(eigenValues[idx]), 10)
    eigenVectors = np.around(
        eigenVectors[:, idx] / eigenVectors[0, idx], 10)
    return (eigenValues, eigenVectors)


def findeigenvec(neigval, eigenvecs):
    '''
        Function for finding particular eigenvector to corresponding eigenvalue's number
    '''

    return eigenvecs[:, neigval]


def modalmatrix(eigval, eigvecs, MassM):
    '''
        Function that creates modal matrix of the system. Inputs are vector of eigenvalues, eigenvectors and mass matrix. Outpu is modal matrix, where each column is normalized eigenvector through modal transform.
    '''

    factors = np.zeros(len(eigval))
    modalM = np.zeros((len(eigvecs[:, 0]), len(eigval)))
    for i in range(len(eigval)):
        factors[i] = np.sqrt(findeigenvec(i, eigvecs).dot(
            MassM.dot(findeigenvec(i, eigvecs))))
        modalM[:, i] = eigvecs[:, i] / factors[i]
    return(np.around(modalM, 10))
