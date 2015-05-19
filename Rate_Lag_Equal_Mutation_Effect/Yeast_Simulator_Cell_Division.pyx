from __future__ import division
import sys
import random
from copy import deepcopy
from heapq import *
from stdlib cimport *
from numpy import *
import numpy as np
cimport numpy as np
cimport cython
import gc
import re  

############################################################################################
# User defined
############################################################################################
import Yeast_Simulator
import Yeast_Simulator_Subprograms
import Yeast_Simulator_Mutations

cdef double mutation_rate = Yeast_Simulator.mutation_rate
cdef double beneficial_mutation_rate = Yeast_Simulator.beneficial_mutation_rate
cdef int n_generations = Yeast_Simulator.number_of_generations
# Cython Imports
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

##############################################################################################
# Each cell divides at differnt speed. 
##############################################################################################
@cython.boundscheck(False)
def asymmetrical_cell_division(genotypes,cell_groups,PM_fitness,mutation_fitness_Lag,
                                PM,n_generations,n_BN,n_individuals,n_mut,point_mutation,text,PM_Track,NN,M_RT,LT,M_LT,CRT,CLT):
    # Cython variable Declarations
    cdef int dps,i,time,tp
    cdef int small
    cdef int k,h,n_cells, P_cdt, C_cdt

    import heapq
    groups = cell_groups.keys()
    groups.sort()
    time = heapq.heappop(groups) 

    # For Fixation analysis
    PM_File = open("PM.txt","w")
    
    # Number of mutations.
    # n_mut = 0
    # desired population size.
    dps = n_individuals * 2**n_generations
    # n_gene_dup = 0
    LL = 0
    #print "Total Population:",dps
    for i in xrange(n_BN):
        #print "Bottleneck :", i
        tp = n_individuals
        # fastest cell division time
        small = genotypes[0][1]
        # slowest cell division time
        big = genotypes[0][1]
        while tp < dps: 
            # Iterating through the cell groups.
            n_cells = len(cell_groups[time])
            for h in xrange(n_cells):
                    k = cell_groups[time][h]
                    # Checking whether the cell can divide or not based on how many times it has divided and Haploids and diploids.
                    #divide = decide_on_division(genotypes[k],div_type)
                    
                    # Change the Lag time to RATE TIME
                    genotypes[k][1] = CRT[k]
                    
                    # Copying the parents attributes to child.
                    genotypes[tp] = genotypes[k]
                    genotypes[tp][2] = 16
                    genotypes[tp][5] = genotypes[k][5]
                    
                    # Update Rate time of the cell
                    CRT[tp] = CRT[k]
                    
                    # Update Lag time of the cell
                    CLT[tp] = CLT[k]
    
                    PM[tp] = PM[k][:]
                    
                    # Introduce mutations
                    if random.random() <= mutation_rate:
                        
                       # Sending mother and daughter cells to introduce mutations. 
                       genotypes[k],genotypes[tp],PM, CRT,CLT = Yeast_Simulator_Mutations.introduce_mutation(genotypes[k],genotypes[tp],
                                                              beneficial_mutation_rate,k,tp,point_mutation,n_mut,PM,PM_fitness,mutation_fitness_Lag,CRT,CLT)
                       n_mut += 1

                    # Reduce the Age for the parent cell
                    genotypes[k][2] = reduce_age(genotypes[k][2])
                    
                    # Updating  the next cell division time based on the cell division time
                    # for both mother and daughter cells.                    
                    (genotypes[k][0],genotypes[tp][0]) = update_next_CDT(genotypes[k][0],genotypes[k][1],genotypes[tp][0],genotypes[tp][1])

                    # Group the cells from current time point
                    P_cdt = genotypes[k][0]; C_cdt = genotypes[tp][0]

                    
                    if P_cdt == C_cdt:               
                       # Grouping the mother cells
                       if cell_groups.has_key(P_cdt):
                          cell_groups[P_cdt].append(k)
                          cell_groups[C_cdt].append(tp)
                       else:
                          cell_groups[P_cdt] = [k]
                          cell_groups[C_cdt].append(tp)
                          heapq.heappush(groups,P_cdt)
                    else:
                       # Grouping the mother cells
                       if cell_groups.has_key(P_cdt):
                          cell_groups[P_cdt].append(k)
                       else:
                          cell_groups[P_cdt] = [k]
                          heapq.heappush(groups,P_cdt)                    
                       # Grouping the daughter cells.
                       if cell_groups.has_key(C_cdt):
                          cell_groups[C_cdt].append(tp)
                       else:
                          cell_groups[C_cdt] = [tp]
                          heapq.heappush(groups,C_cdt)                    
                     
                    tp = cell_count(tp)
                    # Come out of the loop once its reached the desired population size.
                    if tp == dps:
                       break
                    #else:
                    #    genotypes[k][0] += genotypes[k][1]

            del cell_groups[time]
            time = heapq.heappop(groups)

        # Sampling - Bottleneck
        if i != n_BN-1:
            # Random smapling            
            (genotypes,cell_groups,PM,CRT,CLT) = sampling_individuals(genotypes,n_individuals,i,PM,PM_File,n_BN,LT,M_LT,CRT,CLT)
            
        else:
            (genotypes1,cell_groups,PM,CRT,CLT) = sampling_individuals(genotypes,n_individuals,i,PM,PM_File,n_BN,LT,M_LT,CRT,CLT)
            sys.exit()
            
            return genotypes1,PM,PM_fitness,PM_strand,n_mut

        # Its a small mess up. Have to find a better way.
        groups = []
        groups = cell_groups.keys()
        time = heapq.heappop(groups)


#################################################################################################
# Random Sampling the desired number of individuals after certain number of generations.
#################################################################################################
cdef sampling_individuals(genotypes,int n_individuals,int i,PM1,PM_File,n_BN,LT,M_LT,CRT,CLT):

    CRT_Temp = {}
    CLT_Temp = {}
    NBN = i
    # Making the Point mutation and gene duplication Dict of list 
    size = len(genotypes)
    PM = {}; 
    for i1 in range(size):
        PM[i1] = []
    
    # Sampling the number of IDs from genotypes dict.
    cdef int n_cells,j
    cdef float total_time = 0
    size = len(genotypes)
    
    #sampled_ids = random.sample(genotypes,n_individuals)
    sampled_ids = random.random_integers(size-1,size=n_individuals)
    sampled_group = ones((size,6),dtype=int32)
    cell_mutations_temp = []
    cell_gene_duplications_temp = []
    cell_groups = {}
    # Choosing the sampled individuals.
    n_cells = len(sampled_ids)
    n_alive_cells = 0
    for j in xrange(n_cells):
        k = sampled_ids[j]
        if genotypes[k][2] != 0:
            total_time += genotypes[k][1]
            n_alive_cells += 1
            
        #if genotypes[k][5] == 2:
        sampled_group[j] = genotypes[k]
        sampled_group[j][1] = CLT[k]
        sampled_group[j][0] = sampled_group[j][0] + CLT[k]
        CRT_Temp[j] = CRT[k]
        CLT_Temp[j] = CLT[k]
        
        # Copying PM and GD info of cells
        PM[j] = PM1[k][:]
        
        
        # Write it in a file
        #if i == n_BN -1:
        #    print_in_file(PM_File,PM[j],i)
            
        # Storing the mutations in sampled cells.
        # cell_mutations_temp.append(cell_mutations[k])
        # Storing the gene duplications in sampled cells.
        # cell_gene_duplications_temp.append(cell_gene_duplications[k])
        
        # Group the cells based on their cell division time.
        kl = int(genotypes[k][0])
        if cell_groups.has_key(kl):
            cell_groups[kl].append(j)
        else:
            cell_groups[kl] = [j]
            
    Mean_CLT = 0
    Mean_CRT = 0
    for i in xrange(len(CLT)):
        Mean_CLT += CLT[i]
        Mean_CRT += CRT[i]
    
    print NBN+1,"\t",Mean_CLT/len(CLT),"\t","Lag"
    print NBN+1,"\t",Mean_CRT/len(CRT),"\t","Rate"
    
    return sampled_group,cell_groups,PM,CRT_Temp,CLT_Temp

           
###############################################################################################
# To update the next cell division time.
###############################################################################################
cdef inline update_next_CDT(int a1, int a2, int b1, int b2):
     a1 = a1 + a2
     b1 = b1 + b2
     return a1,b1
     
cdef inline reduce_age(int a):
     a = a - 1
     return a
     
cdef inline cell_count(int tp):
     tp = tp + 1
     
     return tp
