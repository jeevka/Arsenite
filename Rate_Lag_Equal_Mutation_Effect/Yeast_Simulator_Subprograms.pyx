####################################################################################################
# This module contains the common subroutines for the Yeast Simulator.
####################################################################################################
####################################################################################################
# Importing Python Modules
####################################################################################################
from __future__ import division
import sys
import random
import operator
from copy import copy
from copy import deepcopy
import math
import re
import numpy as np
from heapq import heapify
from numpy import *
#####################################################################################################
# Importing user defined modules.
#####################################################################################################
import Yeast_Simulator
import Yeast_Simulator_Mutations

################################## Intialization ####################################################
# Intialization of the individuals before the cell division starts.
# 1. Initial population contains cell at different ages. So the cells divide at different speeds.
# 2. 50% cells are virgin cells. 25% are one time divided. 12.5% are two times divided. and etc ...
#####################################################################################################
def initialize_population(n_individuals,n_generations,LT,RT,M_LT):
    cdef int i,n
    size = n_individuals * 2**n_generations
    genotypes = ones((size,6),dtype=int32)
    import random
    population_ratio = {}
    population_ratio[0] = 0
    population_percentage = n_individuals
    #genotypes = {}
    cell_mutations = []
    cell_gene_deletions = []
    cell_gene_duplications = []
    a = []
    
    # Cell RATE TIME
    CRT = {}
    
    # Cell LAG TIME
    CLT = {}
    
    cell_max_age = 16
    # Getting cell division time proportions based on age.
    CDT_proportions = [LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT,LT] 

    # Calculating the ratio/distribution of cells with different ages.
    for i in range(1,18):
        population_ratio[i] = population_ratio[i-1] + population_percentage/2**i

    for n in xrange(n_individuals):
        # Initializing the cell division time and number of generations.
        # Choosing between virgin cells and other cells. Because ~50% cells are virgin cells in the begining.
        # Here n_individuals * 100000 is to get more accurate number upto 5 precision so that the ratio will be better even with small numbers.
        
        ran_max = n_individuals * 100000
        n_ran = (random.randrange(1,ran_max+1))/100000
        n_generations = 16
        
        genotypes[n] = [LT,LT,cell_max_age,0,0,1]
        
        CRT[n] = RT
        CLT[n] = LT
        
    NN = random.randrange(n_individuals)
    #genotypes[NN][1] = M_LT
    #genotypes[NN][0] = M_LT 
    #genotypes[NN][5] = 2
    
    # Making the Point mutation and gene duplication Dict of list 
    PM = {}; 
    for i in range(size):
        PM[i] = []
        
        
    return genotypes,PM,NN,CRT,CLT

#############################################################################################
# Group the cells based on their cell division time.
#############################################################################################
def group_cells(genotypes,n_individuals):
    cell_group = {}
    n_individuals = int(n_individuals)
    for i in xrange(n_individuals):
        if cell_group.has_key(genotypes[i][0]):
            cell_group[genotypes[i][0]].append(i) 
        else:
            cell_group[genotypes[i][0]] = [i]
            
    return cell_group
