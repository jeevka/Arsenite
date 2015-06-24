"""

Example:
# 1. A map distance of 1 Morgan gives a recombination fraction of 0.43
>>> G = array([[0,1],[0,1]])    # simple genotype matrix for two loci
>>> rate = invHaldane(1)        # inverse mapping function of 1 Morgan 
>>> N = 1000                   
>>> nonparental = 0
>>> random.seed(seed=0)
>>> for i in range(N): nonparental += min(sum(recombine(G,rate)))
>>> float(nonparental)/N                        # recombination fraction
0.434
"""

from scipy import mod, cumsum, log, exp
from numpy import random, array, vstack, column_stack, ones, zeros

def Haldane(theta):
    """
    Haldane's mapping function for converting recombination
    fractions to map distances, assuming no interference.

    input:  theta - recombination fraction
    output: map distance (Morgan)
 
    Example:
    >>> Haldane(0.1)
    0.11157177565710485
    >>> Haldane(0.4)
    0.80471895621705025
    """
    return -0.5*log(1-2*theta)

def invHaldane(d):
    """
    The inverse of Haldane's mapping function for converting map distances
    recombination fractions to, assuming no interference.

    input:  d - map distance (Morgan)
    output: recombination fraction
   
    Example:
    >>> invHaldane(0.1)
    0.09063462346100909
    >>> invHaldane(0.4)
    0.27533551794138922
    """
    return 0.5*(1-exp(-2*d))

def recombine(G,rates):
    """
    Performs recombination of genotype matrix G corresponding to given rates

    input:      G       Nx2 matrix of integers where each columns give
                        haplotypes for N markers
                rates   (N-1)x1 vector of floats with recombination rates

    Element rates[i] is the rate of recombination between markers i and i+1

    Example:
    #>>> from numpy import random, ones, zeros, column_stack
    >>> G = column_stack((zeros((50,1)),ones((50,1))))
    >>> rates = 0.2*ones((G.shape[0],1))
    >>> random.seed(seed=0)
    >>> recombine(G,rates).T
    array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
             0.,  0.,  1.,  0.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,
             0.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0.,
             0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.,  0.,  0.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
             1.,  1.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,
             1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.,
             1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0.,  1.,  1.]])
   
    """
   
    #draw uniform numbers for all marker interval
    crossover = mod(cumsum(random.uniform(size=(G.shape[0],1))<rates),2)
   
    #generate recombined G
    recG = array([ [g[c],g[c-1]] for g,c in zip(G[1:G.shape[0],:],crossover) ])
    recG = vstack((G[0,:],recG))
    
    return recG
 
if __name__ == "__main__":
    import doctest
    doctest.testmod()
