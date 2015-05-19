import sys
import random
import numpy as np
import scipy.stats

import Yeast_Simulator_Subprograms

#######################################################################
# This program is the central command of the whole process
#######################################################################
def create_cells(n_individuals,min_n_PM,max_n_PM,min_n_GD,max_n_GD):
    #PM_shape = 0.25; PM_scale = 1
    PM_shape = 0.5; PM_scale = 1
    #GD_shape = 0.5; GD_scale = 1
    GD_shape = 1.5; GD_scale = 1
    
    # Creating the dummy cells.
    genotypes,n_mut,n_gene_dup = create_dummy_cells(n_individuals,min_n_PM,max_n_PM,min_n_GD,max_n_GD)
    
    # Create the haplotypes based on Mut frequencies from sequence data
    create_haplotypes("Point_Mutations.txt",max_n_PM)
    #print "###########################"
    #print set1_min,set1_max
    #print set2_min,set2_max
    #print set3_min,set3_max
    #print set4_min,set4_max
    #print "###########################"
    ### Point mutations ###
    # Point mutation profiles for cells.
    PM_fitness = create_genetic_profiles(genotypes,set4_max+1,PM_shape,PM_scale)
    #print len(PM_fitness)
    #sys.exit()
    # Assign mutations to cells. 3 is the subscript for the point mutations.
    PM = assign_mutations_to_cells(genotypes,PM_fitness,n_individuals,3)
    
    ### Gene duplications ###
    # Create the haplotypes based on Mut frequencies from sequence data
    create_haplotypes("Gene_Duplications.txt",max_n_GD)
    #print "###########################"
    #print set1_min,set1_max
    #print set2_min,set2_max
    #print set3_min,set3_max
    #print set4_min,set4_max
    #print "###########################"
    # Gene duplication profiles for cells.
    GD_fitness = create_genetic_profiles(genotypes,set4_max+1,GD_shape,GD_scale)
    
    # Assign gene duplications to cells. 4 is the subscript for the gene duplication
    GD = assign_mutations_to_cells(genotypes,GD_fitness,n_individuals,4)
    
    # Calculate Cell division time from Genetic profile
    genotypes = calculate_cell_division_time(genotypes,PM,PM_fitness,GD,GD_fitness,n_individuals)
    
    #print "###############################################"
    #print "Thank You Lord"
    #print "###############################################"
    
    return genotypes,PM,PM_fitness,GD,GD_fitness,max_n_PM,max_n_GD
    
#######################################################################
# Create dummy cells with basic properties
#######################################################################
def create_dummy_cells(n_individuals,min_n_PM,max_n_PM,min_n_GD,max_n_GD):
    dummy_array = np.array([180,180,16,0,0])
    genotypes = np.array(np.resize(dummy_array,(n_individuals,5)))
    tot_n_mut = 0; tot_n_GD = 0
    
    # Decising the number of mutations and gene duplications
    for i in xrange(n_individuals):        
        n_PM = random.randrange(min_n_PM,max_n_PM)
        n_GD = random.randrange(min_n_GD,max_n_GD)
        genotypes[i][3] = n_PM
        genotypes[i][4] = n_GD
        tot_n_mut = tot_n_mut + n_PM
        tot_n_GD = tot_n_GD + n_GD
    
    return genotypes,tot_n_mut,tot_n_GD

#######################################################################
# Create dummy cells with basic properties
#######################################################################
def create_genetic_profiles(genotypes,n,shape,scale):
    
    # Create the haplotype for each mutation
    # assign_haplotype(n,file_name)
    
    # Assign Fitness to the mutations
    fitness = assign_fitness(n,shape,scale)
    
    # Assign fitness type
    # fitness = assign_fitness_type(fitness)
    
    return fitness

###################################################################################################
# Create the Haplotypes for the cells based on the mutation freqency from seq data
###################################################################################################
def create_haplotypes(file_name, mut_tot):
    # Needs to be Global for later use
    global set1_min, set1_max, set2_min, set2_max
    global set3_min, set3_max, set4_min, set4_max
    global pro1, pro2, pro3, pro4
    
    # Frequency of the different sets of mutations
    set1 = 1.0; set2= 0.75; set3 = 0.50; set4 = 0.25
    
    # Proportion of the sets
    pro1 = 0.5; pro2 = 0.15; pro3 = 0.075; pro4 = 0.275
    
    # Calculate the number of mutations in each cell and the sample size to sample for particular set
    n_set1,size_set1_sample = n_mut_for_set(mut_tot,set1,pro1)
    n_set2,size_set2_sample = n_mut_for_set(mut_tot,set2,pro2)
    n_set3,size_set3_sample = n_mut_for_set(mut_tot,set3,pro3)
    n_set4,size_set4_sample = n_mut_for_set(mut_tot,set4,pro4)
    
    total_number_mutations = size_set1_sample + size_set2_sample + size_set3_sample + size_set4_sample
    
    # Calling the funtion to create the haplotype.
    assign_haplotype(int(total_number_mutations),file_name)
    
    set1_min = 0; set1_max = size_set1_sample - 1
    set2_min = set1_max + 1; set2_max = set2_min + size_set2_sample - 1
    set3_min = set2_max + 1; set3_max = set3_min + size_set3_sample - 1
    set4_min = set3_max + 1; set4_max = set4_min + size_set4_sample - 1
    
    return 0
    
###################################################################################################
# calculate the number of mutations per cell from the given set and the sample size to sample.
###################################################################################################
def n_mut_for_set(mut_tot,set,pro):
    # total mutation for the given set
    n1 = mut_tot * pro
    # Sample Size
    n2 = n1 * (1/set)
    
    return n1,n2
    
###################################################################################################
# Assigning Fitness effect type
###################################################################################################
def assign_fitness_type(fitness):
    # Fitness type is assumed to follow gumbel distribution
    neutral_min = -0.5;neutral_max = 0.5
    #for i in range(len(fitness)):
        #if fitness[i] >= -0.5 and fitness[i] <= 0.5
    """
    deleterious = 33;neutral = 67; beneficial = 100
    for i in range(len(fitness)):
        n_ran = random.randrange(100)
        if n_ran <= deleterious:
            fitness[i] = fitness[i] * (-1)
        elif n_ran > deleterious and n_ran <= neutral:
            fitness[i] = 0
    """
    
    return fitness

###################################################################################################
# Assigning Fitness effects
###################################################################################################
def assign_fitness(n,shape,scale):
    
    fitness = np.random.gumbel(shape,scale,n)
    fitness_new = []
    for i in xrange(len(fitness)):
        if fitness[i] >= -1.07 and fitness[i] <= 1.07:
            #fitness[i] = 0
            fitness_new.append(0)
        else:
            fitness[i] = fitness[i]/float(32)
            if fitness[i] <= 0.2585:
                fitness_new.append(fitness[i])
    
    return fitness_new

###################################################################################################
# Assigning mutations and gene duplications to cells randomly
###################################################################################################
def assign_mutations_to_cells(genotypes,PM_fitness,n_cells,AS):
    PM = [];
    total_mutations = len(PM_fitness)
    
    for i in xrange(n_cells):
        n_PM = genotypes[i][AS]
        set1,set2,set3,set4 = decide_set_sizes(n_PM)
    
        PM.append([])
    
        # Selecting mutations from differnet sets
        PM = sampling_the_mutations(PM,i,set1,set1_min,set1_max)
        PM = sampling_the_mutations(PM,i,set2,set2_min,set2_max)
        PM = sampling_the_mutations(PM,i,set3,set3_min,set3_max)
        PM = sampling_the_mutations(PM,i,set4,set4_min,set4_max)
        
    return PM

###################################################################################################
# Sampling the mutation ids to add it to cells
################################################################################################### 
def sampling_the_mutations(PM,i,set_size,set_min,set_max):
    
    set_range = range(int(set_min),int(set_max)+1)
    #print "#########################################################3"
    #print set_range
    #print set_size
    ids = random.sample(set_range,set_size)
    
    for j in ids:
        PM[i].append(j)
    
    #print PM[i]
    #print "#########################################################3"
    
    #if i > 10:
    #    sys.exit()
    
    return PM    
        
###################################################################################################
# decide each set size for each cell
###################################################################################################            
def decide_set_sizes(n):
    set1 = int(round(n * pro1))
    set2 = int(round(n * pro2))
    set3 = int(round(n * pro3))
    set4 = int(round(n * pro4))
    
    # Because of rounding, sometime the total becomes more than the n.
    # To avoid that we have to use this simple operataion to change.
    if (set1 + set2 + set3 + set4) > n:
        ran = random.randrange(1,5)
        if ran == 1:
            set1 = set1 - 1
        elif ran == 2:
            set2 = set2 - 1
        elif ran == 3:
            set3 = set3 - 1
        else:
            set4 = set4 - 1
    
    #print "########################"
    #print set1, set2, set3, set4
    #print "########################"
    
    return set1, set2, set3, set4
    
###################################################################################################
# Assigning mutations and gene duplications to cells randomly
###################################################################################################
def assign_mutations_to_cells_1(genotypes,PM_fitness,n_cells,AS):
    PM = [];
    total_mutations = len(PM_fitness)
    for i in xrange(n_cells):
        n_PM = genotypes[i][AS]
        PM.append([])
        for j in xrange(n_PM):
            id = random.randrange(total_mutations)
            PM[i].append(id)
    
    return PM
            
###################################################################################################
# Assigning the Haplotype structure.
###################################################################################################
def assign_haplotype(n,file_name):
    PM = open(file_name,"w")
    # Chromosome numbers and the length.
    chrom = {
             1:230208,2:813178,3:316617,4:1531918,
             5:576869,6:270148,7:1090946,8:562643,
             9:439885,10:745745,11:666454,12:1078175,
             13:924429,14:784333,15:1091289,16:948062
             }
    haplotype_structure = {}
    for i in xrange(n):
        # Choose a chromosome randomly.
        chr_no = random.randint(1,16)
        # Choose a specific position in chromosome randomly.
        bp_position = random.randint(1,chrom[chr_no])
        haplotype_structure[chr_no] = bp_position
        PM.write("{")
        PM.write(str(chr_no))
        PM.write(":")
        PM.write(" ")
        PM.write(str(bp_position))
        PM.write("}")
        PM.write("\n")
    
    PM.close()
    
    return 0

####################################################################################################
# Calculate Cell division time
####################################################################################################
def calculate_cell_division_time(genotypes,PM,PM_fitness,GD,GD_fitness,n):
    for i in xrange(n):
        fitness = 0
        
        # Adding the fitness effects from PM
        for j in xrange(len(PM[i])):
            id = PM[i][j]
            #print "Length:",len(PM_fitness)
            #print "ID:",id                        
            fitness = fitness + PM_fitness[id]
        # Adding the fitness effects from GD
        for j in xrange(len(GD[i])):
            id = GD[i][j]
            #print "Length:",len(GD_fitness)
            #print "ID:",id
            fitness = fitness + GD_fitness[id]        
        
        # Calculate the new CDT based on fitness effects.
        current_CDT = genotypes[i][0]
        genotypes[i][0] = Yeast_Simulator_Subprograms.calculate_cell_division_time(current_CDT,fitness)
        genotypes[i][1] = genotypes[i][0]
        
        if genotypes[i][0] < 0:
            print "Here I am "
            print current_CDT,fitness
            print genotypes[i]
            sys.exit()
                
    return genotypes

#####################################################################################################
#####################################################################################################
