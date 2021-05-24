from scipy import linalg as la
import numpy as np


def eigen(A, B):
    '''
    Výpočet a seřazení vl. čísel (vl. frekvencí) a vl.tvarů od nejmenšího po největší
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
    Funkce ke hledání vl. tvaru k příslušnému zadanému vl. číslu - když chci znát vl. vektor k 1. vl. tvaru zadám 1 jako neigval a eigenvecs je matice, kde jsou seřazeny do sloupců vl. tvary, např. tak, jak vypisuje funkce eigen
    '''

    return eigenvecs[:, neigval]


def modalmatrix(eigval, eigvecs, MassM):
    '''
    Funkce tvoří modální matici systému. Vstupem jsou: vektor vl. čísel eigval, matice vl. vektorů eigvecs a matice hmotnosti MassM. Výstupem je modální matice, kde každý sloupec je normalizovaný vl. vektor pomocí modální transformace, příslušný dané vl. frekvenci.
    '''

    factors = np.zeros(len(eigval))
    modalM = np.zeros((len(eigvecs[:, 0]), len(eigval)))
    for i in range(len(eigval)):
        factors[i] = np.sqrt(findeigenvec(i, eigvecs).dot(
            MassM.dot(findeigenvec(i, eigvecs))))
        modalM[:, i] = eigvecs[:, i] / factors[i]
    return(np.around(modalM, 10))
