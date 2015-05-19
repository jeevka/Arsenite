###################################################################################################################
########################################### IMPORTANT VARIABLES ###################################################
# Genotypes structure goes like this : {id:[next division time, Cell division time, No. of generations,}
# genotypes - contains all the information about the each individual
# mutations - All the mutation structure, beneficial_mutations - contains the IDs of beneficial mutations
# Cell_mutations - cell id and the mutations it carries, gene_duplications - structures of gene duplications.
# cell_gene_duplications - cell id and the gene duplications, gene_deletions - structures of gene deletions
# cell_gene_deletions - cell id and the gene deletions ids, mutations_fitness - mutation id and fitness
# gene_duplication_fitness - gene duplication ID and fitness, gene_deletion_fitness - gene deletion id and fitness.
####################################################################################################################
####################################################################################################################
from __future__ import division
import sys
import cProfile
import scipy.stats
import random
from copy import deepcopy
from numpy import *
import datetime
################################# User defind modules ######################################################
import Yeast_Simulator_Subprograms
#import Yeast_Simulator_Cell_Division
############################################################################################################
####################### Global Variables - Constants ###############################################################
# Refer the coding documents for all the values use here.
number_of_bottlenecks = 50

# Number of generations before bottleneck.
# Which means the average number of cell divisions
# for the cells in initial population.
number_of_generations = 5

# Task ID is used for running 100s of repetitions.
Task_id = sys.argv[1]
mutation_rate =  float(sys.argv[2])
mutation_shape = float(sys.argv[3])
mutation_scale = float(sys.argv[4])
beneficial_mutation_rate = float(sys.argv[5])
fitness_affecting_mutation_rate = float(sys.argv[6])

# WT LAG TIME
LT = 805

# WT RATE TIME
RT = 163

# Mut LAG TIME
M_LT = 271

# Mut RATE TIME
M_RT = 126

# Number of individuals.
n_individuals = 1000

# To store the wildtype cells.
wildtype_cells = {}

# To store the evolved cells.
evolved_cells = {}


##############################################################################################################
####                         Global variables                                                        #########
##############################################################################################################
# Genotypes is a dict which contains all the details of each cell.
genotypes = {} 

# Contains all the mutations
mutations = []

# Contains beneficial mutation's fitness.
# Index will be the same for beneficial_mutations variable.
mutations_fitness = []

# Cell group based on their cell division time.
cell_groups = {}

# desired population size.
dps = n_individuals * 2**number_of_generations

#############################################################################################################
# Subprograms
#############################################################################################################
def make_fitness_proportions(fitness,MP):
    import random as r
    n1 = int(len(fitness) * (100-MP)/100)
    n2 = len(fitness)
    pop = xrange(n2)
    sampled = r.sample(pop,n1)
    for i in sampled:
        fitness[i] *= -1

    return fitness

def covert_to_minutes_Rate(fitness):
    fitness_mins = []
    for i in fitness:
        fitness_mins.append(calculate_cell_division_time_Rate(i))
    
    return fitness_mins

def covert_to_minutes_Lag(fitness):
    fitness_mins = []
    for i in fitness:
        fitness_mins.append(calculate_cell_division_time_Lag(i))
    
    return fitness_mins

def calculate_cell_division_time_Rate(fitness):
    cell_division_max = 163
    cell_division_min = 126
    
    change_CDT = (fitness/(1+fitness)) *  cell_division_max
            
    return change_CDT

def calculate_cell_division_time_Lag(fitness):
    cell_division_max = 805
    cell_division_min = 271
    
    change_CDT = (fitness/(1+fitness)) *  cell_division_max
        
    return change_CDT

def truncate_fitness_effects(fitness,trunc):
    new_fitness = []
    for i in fitness:
        if i <= trunc:
            new_fitness.append(i)
    
    return new_fitness


def calculate_cell_division_time(current_CDT,alpha):
    # Cell Division time to "Alpha"-refer the coding documents.
    cell_division_max = 163
    cell_division_min = 126
    CDT_change = 0
    if alpha > 0:
        CDT_change = (current_CDT - cell_division_min)/(cell_division_max - cell_division_min)
        CDT_change *= alpha
        new_CDT = current_CDT - CDT_change
    else:
        new_CDT = current_CDT + abs(alpha)

    return  new_CDT

def calculate_Lag_time(current_CDT,alpha):
    # Cell Division time to "Alpha"-refer the coding documents.
    cell_division_max = 805
    cell_division_min = 271
    CDT_change = 0
    if alpha > 0:
        CDT_change = (current_CDT - cell_division_min)/(cell_division_max - cell_division_min)
        CDT_change *= alpha
        new_CDT = current_CDT - CDT_change
    else:
        new_CDT = current_CDT + abs(alpha)
    
    #print current_CDT
    #print new_CDT
    
    return  new_CDT
    
    
#############################################################################################################
#   Initializing the Fitnesss effects and haplotype of genetic variants
####################################    #########################################################################
# Fitness effects of mutations for RATE
mutation_fitness_random_fitness = scipy.stats.gamma.rvs(mutation_shape,loc=0,scale=float(1)/float(mutation_scale),size=dps*0.25* number_of_bottlenecks)
mutation_fitness_random_fitness_Rate = covert_to_minutes_Rate(mutation_fitness_random_fitness)
mutation_fitness_random_fitness_Rate = truncate_fitness_effects(mutation_fitness_random_fitness_Rate,36)
mutation_fitness_Rate  = make_fitness_proportions(mutation_fitness_random_fitness_Rate,beneficial_mutation_rate)

# Fitness effects of mutations for LAG
mutation_fitness_random_fitness = scipy.stats.gamma.rvs(mutation_shape,loc=0,scale=float(1)/float(mutation_scale),size=dps*0.25* number_of_bottlenecks)
mutation_fitness_random_fitness_Lag = covert_to_minutes_Lag(mutation_fitness_random_fitness)
mutation_fitness_random_fitness_Lag = truncate_fitness_effects(mutation_fitness_random_fitness_Lag,533)
mutation_fitness_Lag  = make_fitness_proportions(mutation_fitness_random_fitness_Lag,beneficial_mutation_rate)

#############################################################################################################
# Main Program: Contains selection experiment simulation as well as Automated Backcrossings
#############################################################################################################
def Yeast_lab():
    import Yeast_Simulator_Subprograms

    print "Param:",Task_id,"\t",mutation_rate,"\t",mutation_shape,"\t",mutation_scale,"\t",beneficial_mutation_rate,"\t",fitness_affecting_mutation_rate
    
    print 0,"\t",805,"\t","Lag"
    print 0,"\t",163,"\t","Rate"

    # Initialize the cell Population
    genotypes,PM,NN,CRT,CLT = Yeast_Simulator_Subprograms.initialize_population(n_individuals,number_of_generations,LT,RT,M_LT)    
    
    # Grouping the cells based on their cell division time.
    cell_groups =  Yeast_Simulator_Subprograms.group_cells(genotypes,n_individuals)

    # Local variables for the subroutine.     
    #mutations_fitness = {}; PM_strand = {};
    
    # Opening files for storing the HAPLOTYPE OF POINT MUTATIONS AND GENE DUPLICATIONS.
    point_mutation = open("Point_Mutations.txt","w")
    Fname = "PM_Track_" + str(Task_id) + ".csv"
    PM_Track = open(Fname,"w")
    
    import Yeast_Simulator_Cell_Division
    n_mut = 0; 
    
    #ploidy = 1; div_type = 1; s_size = n_individuals; fix_cal = 1
    
    # Calling the function to allow the cells to divide. Mutations, gene deletions and gene duplications will be introduced during cell division.
    (genotypes,PM,PM_fitness,PM_strand,n_mut) = Yeast_Simulator_Cell_Division.asymmetrical_cell_division(genotypes,cell_groups,
                    mutation_fitness_Rate,mutation_fitness_Lag,PM,number_of_generations,number_of_bottlenecks,n_individuals,n_mut,point_mutation,"ASE",PM_Track,NN,M_RT,LT,M_LT,CRT,CLT)
                    
    print "END_OF_ASE"
    PM_Track.close()
    # Closing the HAPLOTYPE files
    
    sys.exit()
#############################################################################################################
#############################################################################################################

#cProfile.run("Yeast_lab()")


# This is to run this module normally like "python Yeast_Simulator.py"
if __name__ == "__main__":
    import sys
    Yeast_lab()