import numpy as np
import math
from rootpy.math.physics.vector import Vector2


def sphericity_aplanarity(vects):
    """
    Calculate the sphericity and aplanarity of a list of TLorentzVectors
    """

    # sphericity tensor
    S = np.zeros(shape=(3,3))
    for i in xrange(3):
        for j in xrange(3):
            S[i][j] = sum([vect.Vect()[i] * vect.Vect()[j] for vect in vects])

    # normalize
    norm = sum([vect.P()**2 for vect in vects])
    S /= norm
    # determine eigenvalues and eigenvectors
    eigvals, eigvects = np.linalg.eig(S)
    # sort eigenvalues in increasing order
    eigvals = sorted(eigvals)

    # sphericity
    sphericity = (eigvals[0] + eigvals[1]) * 1.5

    # aplanarity
    aplanarity = eigvals[0] * 1.5

    return sphericity, aplanarity


def eta_centrality(eta, etaJ1, etaJ2):
    """
    Calculate the eta centrality score for an object to be between two other objects in eta
    Returns 1 if in dead center
    Returns value smaller than 1/e if object is not between
    """

    center = (etaJ1 + etaJ2) / 2.
    width  = 1. / (etaJ1 - center)**2
    return  math.exp(-width * (eta - center)**2)


def phi_centrality(vec_a, vec_b, vec_c):
    """
    Calculate the phi centrality score for an object to be between two other objects in phi
    Returns sqrt(2) if in dead center
    Returns smaller than 1 if an object is not between
    vec_a and vec_b are the bounds, vec_c is the vector to be tested
    """

    aPhi = vec_a.Phi()
    bPhi = vec_b.Phi()
    cPhi = vec_c.Phi()

    A = math.sin(cPhi - aPhi)/math.sin(bPhi - aPhi)
    B = math.sin(bPhi - cPhi)/math.sin(bPhi - aPhi)

    return (A+B)/math.sqrt(A**2 + B**2)
