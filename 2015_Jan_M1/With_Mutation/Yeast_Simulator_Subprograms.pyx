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

#####################################################################################################
# Global Variables - Importing from main program.
#####################################################################################################
cdef int cell_division_min =  Yeast_Simulator.cell_division_min 
cdef int cell_division_max =  Yeast_Simulator.cell_division_max
cdef int n_individuals = Yeast_Simulator.n_individuals
cdef int number_of_bottlenecks = Yeast_Simulator.number_of_bottlenecks
cdef int dps = Yeast_Simulator.dps
mutation_rate = Yeast_Simulator.mutation_rate
beneficial_mutation_rate = Yeast_Simulator.beneficial_mutation_rate
gene_deletion_rate = Yeast_Simulator.gene_deletion_rate
gene_duplication_rate = Yeast_Simulator.gene_duplication_rate
cdef int cell_max_age = Yeast_Simulator.cell_max_age
fcp = Yeast_Simulator.fcp
mcp = Yeast_Simulator.mcp
scp = Yeast_Simulator.scp
#cell_division_time_proportion_distribution
cdt_pd = Yeast_Simulator.cell_division_time_increase_proportions
decay_constant = Yeast_Simulator.decay_constant
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
    
    # Getting cell division time proportions based on age.
    CDT_proportions  = age_cell_divison_time(cell_division_max)
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
        if n_ran >= 0 and n_ran <= population_ratio[1]:
            genotypes[n] = [LT,LT,cell_max_age,0,0,0]
        elif n_ran > population_ratio[1] and n_ran <= population_ratio[2]:    
            genotypes[n] = [CDT_proportions[0], CDT_proportions[0],cell_max_age-1,0,0,0]
        elif n_ran > population_ratio[2] and n_ran <= population_ratio[3]:
            genotypes[n] = [CDT_proportions[1], CDT_proportions[1],cell_max_age-2,0,0,0]
        elif n_ran > population_ratio[3] and n_ran <= population_ratio[4]:
            genotypes[n] = [CDT_proportions[2],CDT_proportions[2],cell_max_age-3,0,0,0]
        elif n_ran > population_ratio[4] and n_ran <= population_ratio[5]:
            genotypes[n] = [CDT_proportions[3],CDT_proportions[3],cell_max_age-4,0,0,0]
        elif n_ran > population_ratio[5] and n_ran <= population_ratio[6]:
            genotypes[n] = [CDT_proportions[4],CDT_proportions[4],cell_max_age-5,0,0,0]
        elif n_ran > population_ratio[6] and n_ran <= population_ratio[7]:
            genotypes[n] = [CDT_proportions[5],CDT_proportions[5],cell_max_age-6,0,0,0]
        elif n_ran > population_ratio[7] and n_ran <= population_ratio[8]:
            genotypes[n] = [CDT_proportions[6],CDT_proportions[6],cell_max_age-7,0,0,0]
        elif n_ran > population_ratio[8] and n_ran <= population_ratio[9]:
            genotypes[n] = [CDT_proportions[7],CDT_proportions[7],cell_max_age-8,0,0,0]
        elif n_ran > population_ratio[9] and n_ran <= population_ratio[10]:
            genotypes[n] = [CDT_proportions[8],CDT_proportions[8],cell_max_age-9,0,0,0]
        elif n_ran > population_ratio[10] and n_ran <= population_ratio[11]:
            genotypes[n] = [CDT_proportions[9],CDT_proportions[9],cell_max_age-10,0,0,0]
        elif n_ran > population_ratio[11] and n_ran <= population_ratio[13]:
            genotypes[n] = [CDT_proportions[10],CDT_proportions[10],cell_max_age-11,0,0,0]
        elif n_ran > population_ratio[12] and n_ran <= population_ratio[14]:
            genotypes[n] = [CDT_proportions[11],CDT_proportions[11],cell_max_age-12,0,0,0]
        elif n_ran > population_ratio[13] and n_ran <= population_ratio[15]:
            genotypes[n] = [CDT_proportions[12],CDT_proportions[12],cell_max_age-13,0,0,0]
        elif n_ran > population_ratio[14] and n_ran <= population_ratio[16]:
            genotypes[n] = [CDT_proportions[13],CDT_proportions[13],cell_max_age-14,0,0,0]
        elif n_ran > population_ratio[15] and n_ran <= population_ratio[17]:
            genotypes[n] = [CDT_proportions[14],CDT_proportions[14],cell_max_age-15,0,0,0]
        else:
            genotypes[n] = [CDT_proportions[15],CDT_proportions[15],cell_max_age-16,0,0,0]
        
        CRT[n] = 162
        CLT[n] = 805
        
    NN = random.randrange(n_individuals)
    #genotypes[NN][1] = M_LT
    #genotypes[NN][0] = M_LT 
    #genotypes[NN][5] = 2
    
    # Making the Point mutation and gene duplication Dict of list 
    PM = {}; 
    for i in range(size):
        PM[i] = []
            
    return genotypes,PM,NN,CRT,CLT

##############################################################################################
# In Reality 10% of the cells wont survive or wont reproduce. This function will do the killings.
##############################################################################################
def inactivate_cells(genotypes,percentage):
    inactive_ids = {}
    n_inactive_cells = int(len(genotypes) * (percentage/100))
    inactive_ids = random.sample(genotypes,n_inactive_cells)    
    for i in range(0,len(inactive_ids)):
        id = inactive_ids[i]
        genotypes[id][2] = 0
    return genotypes

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

#################################################################################################
# Updating the next cell divison time and age
#################################################################################################
def update_CDT_age(parent,daughter):
    # Updating the next cell division time
    parent[0] = int(parent[0]) + int(parent[1])
    daughter[0] = int(daughter[0]) + int(daughter[1])
    # Updating the age
    parent[2] = parent[2] - 1 
    daughter[2] = daughter[2] - 1
    return parent,daughter

#################################################################################################
# Random Sampling the desired number of individuals after certain number of generations.
#################################################################################################
def sampling_individuals(genotypes,n_individuals,cell_mutations,cell_gene_duplications):
    # Sampling the number of IDs from genotypes dict.
    sampled_ids = random.sample(genotypes,n_individuals)        
    sampled_group = {}
    cell_mutations_temp = []
    #cell_gene_deletions_temp = []
    cell_gene_duplications_temp = []
    cell_groups = {}

    # Choosing the sampled individuals.
    for j in xrange(len(sampled_ids)):
        k = sampled_ids[j]
       
        sampled_group[j] = genotypes[k][:]
        # Storing the mutations in sampled cells.
        cell_mutations_temp.append(cell_mutations[k])
        # Storing the gene deletions in sampled cells.
        #cell_gene_deletions_temp.append(cell_gene_deletions[k])
        # Storing the gene duplications in sampled cells.
        cell_gene_duplications_temp.append(cell_gene_duplications[k])
        # Group the cells based on their cell division time.
        kl = int(genotypes[k][0])
        #print "KLKL:",kl
        if cell_groups.has_key(kl):
            cell_groups[kl].append(j)
        else:
            cell_groups[kl] = [j]

    #print "Sampling :",len(cell_mutations)
    return sampled_group,cell_groups,cell_mutations_temp,cell_gene_duplications_temp

#################################################################################################
# non-Random : Sampling the desired number of individuals after certain number of generations.
#################################################################################################
def sampling_individuals_test(genotypes,n_individuals,cell_mutations,cell_gene_duplications,small,big):
    # Sampling the number of IDs from genotypes dict.
    sampled_ids = group_the_cells(genotypes,n_individuals,small,big)

    sampled_group = {}
    cell_mutations_temp = []
    #cell_gene_deletions_temp = []
    cell_gene_duplications_temp = []
    cell_groups = {}

    # Choosing the sampled individuals.
    for j in xrange(len(sampled_ids)):
        k = sampled_ids[j]
        sampled_group[j] = genotypes[k][:]
        # Storing the mutations in sampled cells.
        cell_mutations_temp.append(cell_mutations[k])
        # Storing the gene deletions in sampled cells.
        #cell_gene_deletions_temp.append(cell_gene_deletions[k])
        # Storing the gene duplications in sampled cells.
        cell_gene_duplications_temp.append(cell_gene_duplications[k])
        # Group the cells based on their cell division time.
        kl = int(genotypes[k][0])
        #print "KLKL:",kl
        if cell_groups.has_key(kl):
            cell_groups[kl].append(j)
        else:
            cell_groups[kl] = [j]
    #print "Sampling :",len(cell_groups)
    return sampled_group,cell_groups,cell_mutations_temp,cell_gene_duplications_temp


##############################################################################################
# This function will recalculate the cell division time based on its Age
##############################################################################################
def age_based_division_time_change(mother,daughter):
    age = 16 - mother[2]
    mother[1] = mother[1] + (cdt_pd[age-1] * mother[1])

    return mother,daughter
    
##############################################################################################
# This function will break the whole population into 3 major groups based on their division time.
##############################################################################################
def group_the_cells(genotypes,n,small,big):
    #fa = 40;me = 20;sl = 40
    break_point = (big - small)/3
    p1 = small + break_point
    p2 = p1 + break_point
    p3 = p2 + break_point
    t = []
    # Faster dividing cells
    faster = {}
    # Medium dividing cells
    medium = {}
    # Slower dividing cells
    slower = {}
    for i in xrange(0,len(genotypes)):
        if genotypes[i][1] < p1:
            faster[i] = genotypes[i][1]
        elif genotypes[i][1] < p2:
            medium[i] = genotypes[i][1]
        else:     
            slower[i] = genotypes[i][1]
    # Decide number of cells in each group based on available population size in each group
    if (float(n)*(fcp/100)) > len(faster) or (float(n)*(mcp/100)) > len(medium) or (float(n)*(scp/100)) > len(slower):
        (n1,n2,n3) = resize_the_cell_groups(n,fcp,mcp,scp,faster,medium,slower)
    else:
        n1 = int(float(n) * (float(fcp)/float(100)))
        n2 = int(float(n) * (float(mcp)/float(100)))
        n3 = int(float(n) * (float(scp)/float(100)))
    #print "Size:", len(faster),len(medium),len(slower)
    # Sampling each group seperately
    #print "NNN:",n
    sampled_ids1 = random.sample(faster,n1)
    sampled_ids2 = random.sample(medium,n2)
    sampled_ids3 = random.sample(slower,n3)    
    
    # merging the arrays
    for i in xrange(0,len(sampled_ids2)):
        sampled_ids1.append(sampled_ids2[i])
    for i in xrange(0,len(sampled_ids3)):
        sampled_ids1.append(sampled_ids3[i])
    
    return sampled_ids1

###############################################################################################    
# This is to re assign the number of cells to be sampled from each group.
###############################################################################################
def resize_the_cell_groups(n,fa,me,sl,faster,medium,slower):
    n1 = int(float(n) * (float(fa)/float(100)))
    n2 = int(float(n) * (float(me)/float(100)))
    n3 = int(float(n) * (float(sl)/float(100)))
    #print "Before:",n1,n2,n3
    #print len(faster),len(medium),len(slower)
    if n1 > len(faster):
       d = n1 - len(faster)
       n1 = n1 - d
       n2 = n2 + d
    if n2 > len(medium):
       d = n2 - len(medium)
       n2 = n2 - d
       # Decide which cells to add
       if n1 > n3:
        n1 = n1 +d
       else:
        n3 = n3 + d
    if n3 > len(slower):
       d = n3 - len(slower)
       n3 = n3 - d
       n2 = n2 + d
    # print "After:",n1,n2,n3   
    return n1,n2,n3

###############################################################################################    
# This is to calculate the mean cell divison time of the cells at each bottlneck.
###############################################################################################
def calculate_mean_cell_division_time(genotypes,output_file,i):
    # "i" is the bottlneck number.
    total_time = 0
    n_cells = 0
    for l in xrange(0,len(genotypes)):
        # Ignore the dead cells.
        if genotypes[l][2] != 0:
            total_time = total_time + genotypes[l][1]
            n_cells = n_cells + 1
    mean_division_time = total_time/n_cells
    print i,mean_division_time
    
    return mean_division_time

##############################################################################################
# This is to calculate the Cell Division time based on the model we developed to do.
# Look into the coding documents.
##############################################################################################
def calculate_cell_division_time(cell_division_time,alpha):

    # Cell Division time to "Alpha"-refer the coding documents.
    pro = 1 - (cell_division_max - cell_division_time)/(cell_division_max - cell_division_min)
    alpha_current = math.log(pro)/-(math.log(decay_constant))
    CDT = cell_division_time
    
    # Add the old alpha with new Alpha to calculate the new cell division time.
    if alpha > 0:
        alpha_new = alpha + alpha_current
        beta = cell_division_max - cell_division_min
        natural_log = -(math.log(decay_constant) * (alpha_new))
        CDT = math.ceil(beta * math.exp(natural_log) + cell_division_min)
    elif alpha == 0:
        alpha = 0
    else:
        alpha_new = (alpha + alpha_current)
        beta = cell_division_max - cell_division_min
        natural_log = -(math.log(decay_constant) * (alpha_new))
        CDT = math.ceil(beta * math.exp(natural_log) + cell_division_min)

    return CDT

###############################################################################################
# To find mean number of mutations, gene deletions, gene duplications.
###############################################################################################
def calculate_mean_number_of_genetic_variations(cell_mutations,cell_gene_deletions,cell_gene_duplications):
   
    mean_mutation = 0
    mean_gene_deletion = 0
    mean_gene_duplication = 0
    
    for i in xrange(len(cell_mutations)):
        mean_mutation = mean_mutation + len(cell_mutations[i])
        
    print "Mean Mutation:",(mean_mutation)/float(len(cell_mutations))
    

    for i in xrange(len(cell_gene_duplications)):
        mean_gene_duplication = mean_gene_duplication + len(cell_gene_duplications[i])

    print "Mean Gene Duplication :",(mean_gene_duplication)/float(len(cell_gene_duplications))

    return 0

###############################################################################################
# Calculate the Cell division time increase based on age.
###############################################################################################
def age_cell_divison_time(cell_division_max):
    m = 0.854e12; b = 0.765; r = 4.925e-14
    f = 0; S = 0; g = 0.0421; h = 3.7676e-14
    time = []
    CDT_proportions = []
    # Calculating the cell division time.
    for n in range(0,17):
        CDT = ((m/g) * ((h + (g * S)) * math.exp(g * n)) - h) + b
        time.append(CDT)
        
    # calculating the Percentage of increase.
    for i in range(1,len(time)):
        t1 = time[i] - time[i-1]
        t2 = t1/time[i-1]
        t3 = cell_division_max + int(cell_division_max * t2)
        CDT_proportions.append(t3)
    
    return CDT_proportions
    
###################################################################################################
# Assigning the Haplotype structure.
###################################################################################################
def assign_haplotype(n):
        
    # Chromosome numbers and the length.
    chrom = {
             1:230208,2:813178,3:316617,4:1531918,
             5:576869,6:270148,7:1090946,8:562643,
             9:439885,10:745745,11:666454,12:1078175,
             13:924429,14:784333,15:1091289,16:948062
             }

    haplotype_structure = []
    
    for i in range(n):

        # Choose a chromosome randomly.
        chr_no = random.randint(1,16)

        # Choose a specific position in chromosome randomly.
        bp_position = random.randint(1,chrom[chr_no])
        haplotype_structure.append(str(chr_no) + str(":") + str(bp_position))
        
    return haplotype_structure
#################################################################################################### 