from __future__ import division
import random
import sys
import math
import numpy
import scipy

# Cython Imports
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

###################################################################################################
# Import user defined Modules.
###################################################################################################
import Yeast_Simulator
import Yeast_Simulator_Subprograms
import Yeast_Simulator_Cell_Division
###################################################################################################
# Global variables - Importing from main Modules.
###################################################################################################
fitness_affecting_mutation_rate = Yeast_Simulator.fitness_affecting_mutation_rate

###################################################################################################
# Introducing Mutations Based on Yeast Mutation rate.
###################################################################################################
def introduce_mutation(mother_cell,daughter_cell,beneficial_mutation_rate,k,tp,point_mutation,n_mut,PM,PM_fitness,PM_fitness_Lag,CRT,CLT):

    # Decide whether its fitness altering mutation or neutral.
    # Beneficial Mutation rate is 575 out of 10000 fitness altering mutations.
    # Its 0.088 or 88/1000. Plz refer the coding value documents.
    # From Coding documents: 1.53% of all mutations are Fitness altering mutations.
    
    fitness_rate = 0
    fitness_lag = 0
    # Take Fitness effect for RATE
    if random.randrange(100) < fitness_affecting_mutation_rate:        
        fitness_rate = PM_fitness[n_mut]
        fitness_lag = 0 #PM_fitness_Lag[n_mut]
    else:
        fitness_rate = 0 
    
    """
    # Take Fitness effect for LAG
    if random.randrange(100) < fitness_affecting_mutation_rate:        
        fitness_lag = PM_fitness_Lag[n_mut]
    else:
        fitness_lag = 0
    """
            
    # Choose the cell:Parent or child to have mutation.
    ran_n = c_libc_random() % 0.99 + 0
    if ran_n <= 0.5:
        
        if fitness_rate != 0:
            # Calculate the RATE CDT
            mother_cell[1] = Yeast_Simulator.calculate_cell_division_time(mother_cell[1],fitness_rate)
            CRT[k] = mother_cell[1]
            
        if fitness_lag != 0:
            # Calculate the Lag Time  
            CLT[k] = Yeast_Simulator.calculate_Lag_time(CLT[k],fitness_lag)
            #CLT[k] = mother_cell[1]
            
            mother_cell[5] = 2
            
        point_mutation.write("\n")
        mother_cell[3] += 1
        
        # Appending the mutation number
        PM[k].append(n_mut)
        
    else:        
        if fitness_rate != 0:
            # Calculate the RATE CDT
            daughter_cell[1] = Yeast_Simulator.calculate_cell_division_time(daughter_cell[1],fitness_rate)
            CRT[tp] = daughter_cell[1]
        
        if fitness_lag != 0:
            # Calculate the Lag Time  
            CLT[tp] = Yeast_Simulator.calculate_Lag_time(CLT[tp],fitness_lag)
            #CLT[tp] = daughter_cell[1]        
            
            daughter_cell[5] = 2
        
        point_mutation.write("\n")
        daughter_cell[3] +=  1
       
        # Appending the mutation number
        PM[tp].append(n_mut)
    
    return mother_cell,daughter_cell,PM, CRT,CLT

####################################################################################################    

"""
Reference:
Sarah B et al.  Spontaneous Mutations in Diploid Saccharomyces Cerevisiae: More Beneficial Than Expected. 2004 Genetics.
"""
####################################################################################################    