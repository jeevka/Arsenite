###################################################################################################################
############################################# MAIN PROGRAM ########################################################
###################################################################################################################
# MAIN CONCEPT: Initial population contains n(Ex.n=10) number of yeast cells and each cell divides at different
# time point (Ex. bw 78-108 mins). When the time goes on, each cell divides at its particular time point.
# Some cells grows faster(78 mins) and some are slower(102 mins). Cell growth will be stopped when the desird
# population size is reached. And then n number of cells will be randomly sampled(bottleneck) and let it grow again.
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
import BackCrossing
#import Artificial_Cells
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
# Min and max cell division time for the cells
# cell_division_min = 78
# cell_division_max = 102
# Use this value for stressor medium.
#cell_division_min = 78
# Mean cell division time of the WT cells in normal medium
cell_division_min = 90
# Maximum cell divison time of the cells in Stressor medium.
cell_division_max = 180
# Maximum Number of times a cell can divide.
cell_max_age = 16
# Number of chromosomes.
n_chr = 16
# Number of Genes in Yeast Genome.
n_genes = 5802
# Mutation rate is 4 out of 1000.
Task_id = sys.argv[1]

decay_constant = 2 #float(sys.argv[2])
mutation_rate =  float(sys.argv[2])
mutation_shape = float(sys.argv[3])
mutation_scale = float(sys.argv[4])
# Beneficial Mutation rate is 88 out of 1000 in Genome wide mutation rate.
beneficial_mutation_rate = float(sys.argv[5]) #60 #0.088
fitness_affecting_mutation_rate = float(sys.argv[6])
beneficial_mutation_fitness_min = 5  # 5%
beneficial_mutation_fitness_max = 10 # 10%
# This is rate is per gene per generation.
# gene_duplication_rate = 0.0000034
# gene_duplication_rate = 0.0197268 #float(sys.argv[4])
gene_duplication_rate = 0 #0.0195942 #float(sys.argv[4])
gene_duplication_shape = 0.4 #2 #float(sys.argv[5])
gene_duplication_scale = 0.2 #1 #float(sys.argv[6])
gene_duplication_fitness = 0.5
gene_duplication_beneficial = 100 #float(sys.argv[7])
gene_duplication_neutral = 0 #float(sys.argv[8])
gene_duplication_deleterious = 0 #float(sys.argv[9])
gene_duplication_fitness_truncation = 16 # float(sys.argv[10])
# Non-Random Sasmpling
# Faster cells proportion 
fcp = 40
# Medium cells proportion
mcp = 20
# Slower cells proportion
scp = 40
# Gene Duplication Rate is 1218 out of 100000
# Value is converted from 1 gene rate to 5802 genes.
gene_deletion_rate = 0.0121842 # Out of 10000000
# Number of gene duplications
n_gene_duplications = 0
# Number of gene deletions.
n_gene_deletions = 0
# Beneficial Mutation rate is 575 out of 10000.
# beneficial_mutation_rate = 0.0575
# Number of individuals.
n_individuals = 1000
# cell % of wont viable
cell_wont_survive = 10
# Variables which are used for Testing purpose.
number_of_mitosis_events = 0
total_number_of_mitosis_events = 0
total_population_size = 0

# Number of cells for mitosis is between 20-25
# Number of evolved cells to backcross.
import random
min_sample_size = 20
max_sample_size = 25
n_evolved_cells = random.randrange(min_sample_size,max_sample_size)
# Number of wildtype to backcross.
n_wildtype = n_evolved_cells
n_backcrossing = 10
# Mating Parameters
p_m = 0.90

# To store the wildtype cells.
wildtype_cells = {}
# To store the evolved cells.
evolved_cells = {}

# n_mut will tell the number of mutation events happened.
# tp: total population. tp will tell the population size.

# Age and Cell division time distribution.
# Ex. cell_division_time{age:cellDivisionTime}
"""
cell_division_time = {0:102,1:102,2:98,3:94,
                      4:88,5:84,6:78,7:80,8:82,
                      9:84,10:86,11:88,12:90,
                      13:92,14:94,15:96,16:102}
cell_division_time_distribution = [78,80,82,84,86,
                                   88,90,92,9,94,96,102]
"""
cell_division_time_increase_proportions = [0.0219414690133,0.0223935772857,0.0228448947257,0.0232950247887,0.0237435751105,0.0241901588628,0.0246343960682,0.0250759148643,
                                           0.0255143527094,0.0259493575199,0.0263805887355,0.0268077183016,0.0272304315679,0.0276484280956,0.0280614223718,0.0284691444277]

##############################################################################################################
####                         Global variables                                                        #########
##############################################################################################################
# Genotypes is a dict which contains all the details of each cell.
genotypes = {} 
# Contains all the mutations
mutations = []
# Contains beneficial Mutations
beneficial_mutations = []
# Contains beneficial mutation's fitness.
# Index will be the same for beneficial_mutations variable.
mutations_fitness = []
# Cell ID and the beneficial Mutation IDs
cell_beneficial_mutations = []
# Contains the haplotype of gene duplications
gene_duplications = []
# Contains the gene duplication fitness.
# Index will be the same as gene_duplications variable.
gene_duplications_fitness = []
# Contains the cell ID and gene duplication ids.
cell_gene_duplications = []
# Contains Gene deletion haplotypes
gene_deletions = []
# Contains gene deletion fitness
# Index is same as gene_deletions variable.
gene_deletions_fitness = []
# Contains cell id and the gene deletion IDs
cell_gene_deletions = []
n_individuals_with_mutation = 0
# This variables will be used in Mutation freq calculalation.
bottleneck_positions = {}
bottleneck_positions[0] = 0
# Cell group based on their cell division time.
cell_groups = {}
# Total number of mutation events.
n_mut = 0

cell_mutations = array([])

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
    
def scale_fitness_effects(fitness):
    fitness_scale = make_scales()
    for i in range(len(fitness)):
        id = "%.4f" % round(fitness[i],4)
        fitness[i] = fitness_scale[id]
    
    return fitness
    
def make_scales():
    scale = {}
    alpha = 0
    for i in range(50000):
        P = -math.log(2) * alpha 
        CDT1 = 90 * math.exp(P) + 90
        CDT2 = 180 - 90*alpha
        new_alpha = alpha
        if CDT1 < CDT2:
            new_alpha = find_new_alpha(CDT1,CDT2,alpha)

        alpha1 = "%.4f" % round(alpha,4)
        
        scale[alpha1] = new_alpha
        alpha += 0.0001
        
    return scale


def covert_to_minutes(fitness):
    fitness_mins = []
    for i in fitness:
        fitness_mins.append(calculate_cell_division_time_1(i))
    
    return fitness_mins

def covert_to_minutes_1(fitness):
    fitness_mins = []
    for i in fitness:
        fitness_mins.append(calculate_cell_division_time_2(i))
    
    return fitness_mins

def calculate_cell_division_time_1(alpha):
    cell_division_max = 180
    cell_division_min = 90
    change_CDT = (cell_division_max-cell_division_min) * alpha
    
    if change_CDT >= 90:
        change_CDT = 89
        
    return change_CDT


def calculate_cell_division_time_2(alpha):
    cell_division_max = 660
    cell_division_min = 300
    change_CDT = (cell_division_max-cell_division_min) * alpha
    
    if change_CDT >= 359:
        change_CDT = 359
        
    return change_CDT
    
    
def truncate_fitness_effects(fitness,trunc):
    new_fitness = []
    for i in fitness:
        if i < trunc:
            new_fitness.append(i)
    
    return new_fitness


def calculate_cell_division_time(current_CDT,alpha):
    #print "current_CDT:",current_CDT
    #print "alpha:",alpha
    
    # Cell Division time to "Alpha"-refer the coding documents.
    cell_division_max = 180
    cell_division_min = 126
    CDT_change = 0
    if alpha > 0:
        CDT_change = (current_CDT - cell_division_min)/(cell_division_max - cell_division_min)
        CDT_change *= alpha
        new_CDT = current_CDT - CDT_change
    else:
        new_CDT = current_CDT + abs(alpha)
    
    #print "new_CDT:",new_CDT
    
    return  new_CDT

def calculate_Lag_time(current_CDT,alpha):
    # Cell Division time to "Alpha"-refer the coding documents.
    cell_division_max = 660
    cell_division_min = 300
    CDT_change = 0
    if alpha > 0:
        CDT_change = (current_CDT - cell_division_min)/(cell_division_max - cell_division_min)
        CDT_change *= alpha
        new_CDT = current_CDT - CDT_change
    else:
        new_CDT = current_CDT + abs(alpha)

    return  new_CDT
    
    
#############################################################################################################
#   Initializing the Fitnesss effects and haplotype of genetic variants
#############################################################################################################
mutation_fitness_random_fitness = scipy.stats.gamma.rvs(mutation_shape,loc=0,scale=float(1)/float(mutation_scale),size=dps*0.25* number_of_bottlenecks)
mutation_fitness_random_fitness = scale_fitness_effects(mutation_fitness_random_fitness)
mutation_fitness_random_fitness
mutation_fitness_random_fitness_Rate = covert_to_minutes(mutation_fitness_random_fitness)
mutation_fitness_random_fitness_Lag = covert_to_minutes_1(mutation_fitness_random_fitness)
  
mutation_fitness_random_fitness_Rate = truncate_fitness_effects(mutation_fitness_random_fitness_Rate,54)
mutation_fitness_random_fitness_Lag = truncate_fitness_effects(mutation_fitness_random_fitness_Lag,359)

mutation_fitness_Rate  = make_fitness_proportions(mutation_fitness_random_fitness_Rate,beneficial_mutation_rate)
mutation_fitness_Lag  = make_fitness_proportions(mutation_fitness_random_fitness_Lag,beneficial_mutation_rate)

#############################################################################################################
# Main Program: Contains selection experiment simulation as well as Automated Backcrossings
#############################################################################################################
def Yeast_lab():
    import Yeast_Simulator_Subprograms
    import BackCrossing
    N_HAP1 = 0;N_HAP2 = 0;N_HAP3 = 0

    print "CDT\t",0,"\t",805,"\t","Lag"
    print "CDT\t",0,"\t",162,"\t","Rate"

    # WT LAG TIME
    LT = 805
    
    # WT RATE TIME
    RT = 162
    
    # Mut LAG TIME
    M_LT = 300
    
    # Mut RATE TIME
    M_RT = 126
    
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
    # Cell_Haplotype = {}
    CH =    {}
    N_CH = {}
    MID_Hap = {}
    for i in xrange(n_individuals):
        CH[i] = "W"
        N_CH[i] = 0
    
    # Ploidy = 1 : Haploid, ploidy = 2: Diploid
    ploidy = 1; div_type = 1; s_size = n_individuals; fix_cal = 1
    # Calling the function to allow the cells to divide. Mutations, gene deletions and gene duplications will be introduced during cell division.
    (genotypes,PM,PM_fitness,PM_strand,n_mut) = Yeast_Simulator_Cell_Division.asymmetrical_cell_division(genotypes,cell_groups,
                    mutation_fitness_Rate,mutation_fitness_Lag,PM,number_of_generations,number_of_bottlenecks,n_individuals,n_mut,point_mutation,ploidy,"ASE",div_type,s_size,fix_cal,PM_Track,NN,M_RT,LT,M_LT,CRT,CLT,N_HAP1,N_HAP2,N_HAP3,CH,N_CH,MID_Hap)
            
    print "END_OF_ASE"
    PM_Track.close()
    # Closing the HAPLOTYPE files
    
    sys.exit()
    
    point_mutation.close()
    
    #file_name = "Epistasis_" + str(Task_id) + str(".csv")    
    #epistasis = open(file_name,"w")
    
    # Opening the Point mutation Haplotype files
    f1 = open("Point_Mutations.txt","r")
    PM_hap = f1.readlines()
    f1.close()
    
    # Deciding the number of cells to backcross
    # Extra 10% cells are added since the number of haploids to mate or ~50% not exactly 50%
    Ini_n_evolved_cells = n_individuals
    
    # Number of wildtype to backcross.
    Ini_n_wildtype = Ini_n_evolved_cells
    
    #n_mut = Ini_n_mut; n_gene_dup = Ini_n_gene_dup

    # Selecting Wildtype cells for backcrossing
    # Ini_wildtype_cells = BackCrossing.choose_wild_type_cells(Ini_n_wildtype)
    
    
    # Selecting the highly evolved cells from the population to backcross with founder population.    
    Ini_evolved_cells,Ini_PM_2_BC,Ini_PM_fitness_2_BC,Ini_PM_hap_2_BC = BackCrossing.select_evolved_cells_randomly_to_cross(Ini_n_evolved_cells,genotypes,n_individuals,PM,PM_fitness,PM_hap)    
        
    # Get the allels in the selected populaiton to track their frequency.
    #GDs = BackCrossing.combine_fitness(Ini_GD_2_BC,Ini_GD_fitness_2_BC)
    #PMs = BackCrossing.combine_fitness(Ini_PM_2_BC,Ini_PM_fitness_2_BC)
            
    initial_PMs = Ini_PM_fitness_2_BC.keys()
    backcross_number = 0
    
    # Dummy file
    PMT = open("PM_Track_Dummy.csv","w")
    
    #######################################################################################
    # Repeat the Entire BackCrossing Procedure Number of Times
    #######################################################################################
    for j in range(10):
        print "Repeating the entire BackCrossings:",j
        #print "TIME BEFORE:",j,datetime.datetime.now() 
        
        evolved_cells = deepcopy(Ini_evolved_cells)
        
        PM_2_BC = deepcopy(Ini_PM_2_BC)
        PM_fitness_2_BC = deepcopy(Ini_PM_fitness_2_BC)
        PM_hap_2_BC = deepcopy(Ini_PM_hap_2_BC)        
        # Decide the parameters
        n_BC_generations = 5; N_BN = 1; param_text = "5BN_5GN"
        
        PM_new_ids = []

        n_wildtype = Ini_n_wildtype
        n_evolved_cells = Ini_n_evolved_cells
        
        #######################################################################################################
        # STEP XVI: Repeat Step 1-16 number of times
        #######################################################################################################
        # Number of repeat backcrossing
        for i in xrange(n_backcrossing):
            print "BackCrossing:",i
            
            #######################################################################################################
            # Explantion for some of the indicator variables used in the code
            # fix_cal : 1. Whether fixation calculation is needed or not. 2. Whether PMs and GDs carry fitness effect or not
            # ploidy  : Ploidy nature of the genome. Depends upon ploidy, chr-strand for PMs and GDs are being decided
            # div_type: Whether its, mitotic division for Haploids or Diploids or together
            #######################################################################################################
            
            #######################################################################################################            
            # STEP I and II:  Mix WT HAPLOID and evolved cells(HAPLOID) and Allow 50% of them to mate
            #######################################################################################################
                
            # Back Crossing. Allow the haploid cells to mate and produce diploids,then recombination , and then each diploid cell will divide into four haploids.
            evolved_cells,n_cell,PM_new,PM_fitness_new,PM_hap_new,s_cells,PM_strand = BackCrossing.back_cross(evolved_cells,PM_2_BC,PM_fitness_2_BC,PM_hap_2_BC,initial_PMs,j,i,p_m,PM_strand)
            
            # This FREQ Calculation does it only for the selected cells
            if j == 0 and i == 0:
                #BackCrossing.calculate_frequency_1('GD',"B4BC",j,backcross_number,GD_2_BC,initial_GDs,GD_fitness_2_BC,param_text,GD_new_ids,s_cells)
                BackCrossing.calculate_frequency_1('PM',"B4BC",j,backcross_number,PM_2_BC,initial_PMs,PM_fitness_2_BC,param_text,PM_new_ids,s_cells)
                print "END_OF_B4BC"
            
            #######################################################################################################            
            # STEP III: Allow Mitotic division of Mated and Non-mated cells for "five" generations 
            #######################################################################################################
            
            #LOST1 = BackCrossing.freq_test2(PM_new,initial_PMs,500,PM_fitness_2_BC)
            #print "LOST1:",i,"\t",j,"\t",LOST1
            
            
            # To determine whether its a selection experiment under stress conditions or not.
            # Assumption: During Non-selection experiment, all the cells grow in the same speed.
            genotypes = BackCrossing.prepare_for_non_s_exp(evolved_cells,n_cell)
            
            # Grouping the cells based on their cell division time for "DIFFERENTIAL REPRODUCTION".
            cell_groups =  Yeast_Simulator_Subprograms.group_cells(evolved_cells,n_cell)
            
            # Converting the data format: Changing the EXISTING program may cost lots of time to change.
            # So better to change the format.
            
            # GD_temp = BackCrossing.convert_format(GD_new)
            # GD_fitness_temp = BackCrossing.combine_fitness(GD_new,GD_fitness_new)
            PM_temp = BackCrossing.convert_format(PM_new)
            # PM_fitness_temp = BackCrossing.combine_fitness(PM_new,PM_fitness_new)
            
            # Decide the parameters
            n_BC_generations = 5; N_BN = 1; param_text = "5BN_5GN"
                        
            # Expanding the evolved_cells size
            genotypes_1 = BackCrossing.expand_format(evolved_cells,n_BC_generations,N_BN,n_cell)        
            
            # RE-Opening the HAPLOTYPE files to append the new HAPLOTYPE OF POINT MUTATIONS AND GENE DUPLICATIONS.
            point_mutation = open("Point_Mutations.txt","a")
            
            N_BN = 1; s_size = n_cell
            ploidy = 2; div_type = 1; fix_cal = 0
            #from datetime import datetime            
            # Differential Reproduction for EVOLVED HAPLOIDs, WT HAPLOIDS and DIPLOIDS
            (genotypes_1,PM_2_BC,PM_fitness_2_BC,PM_strand,n_mut) = Yeast_Simulator_Cell_Division.asymmetrical_cell_division(genotypes_1,cell_groups,PM_fitness_new,PM_temp,PM_strand,n_BC_generations,N_BN,n_cell,n_mut,point_mutation,ploidy,"DBC",div_type,s_size,fix_cal,PMT)
            
            # Closing the HAPLOTYPE files
            point_mutation.close()
            
            
            # Opening the Point mutation Haplotype files
            f1 = open("Point_Mutations.txt","r")
            PM_hap_2_BC = f1.readlines()
            f1.close()
            
            #LOST2 = BackCrossing.freq_test3(PM_2_BC,initial_PMs,550,PM_fitness_2_BC)
            #print "LOST2:",i,"\t",j,"\t",LOST2           
            
            #  Opening the Gene Duplication Haplotype files
            # f1 = open("Gene_Duplications.txt","r")
            # GD_hap_2_BC = f1.readlines()
            # f1.close()
            
            # Number of wildtype to backcross.
            n_wildtype = n_individuals/2 
            n_evolved_cells = n_individuals/2
            
            # Checking funtion to fix a bug: its better to keep to check this for future.
            # BackCrossing.check_function(PM_2_BC,PM_fitness_2_BC)
            # BackCrossing.check_function(GD_2_BC,GD_fitness_2_BC)
            # sys.exit()
            #######################################################################################################
            # STEP IV: Subsample 50000 of them and CHOOSE only the DIPLOIDS: THIS STEP CAN BE AVOIDED
            #######################################################################################################            
            
            # Selecting the highly evolved cells from the population to backcross with founder population.
            # evolved_cells,GD_2_BC,GD_fitness_2_BC,GD_hap_2_BC,PM_2_BC,PM_fitness_2_BC,PM_hap_2_BC = BackCrossing.select_evolved_cells_randomly_to_cross(n_evolved_cells,genotypes_1,n_cell,GD_2_BC,GD_fitness_2_BC,GD_hap_2_BC,PM_2_BC,PM_fitness_2_BC,PM_hap_2_BC)

            #######################################################################################################
            # STEP V: Allow the Diploids to grow mitotically for 5 Gs
            #######################################################################################################
    
            # STEP 5.1: Select only the DIPLOID CELLS
            # evolved_cells,GD_2_BC,GD_fitness_2_BC,GD_hap_2_BC,PM_2_BC,PM_fitness_2_BC,PM_hap_2_BC,n_diploids = BackCrossing.select_diploid_cells(evolved_cells,GD_2_BC,GD_fitness_2_BC,GD_hap_2_BC,PM_2_BC,PM_fitness_2_BC,PM_hap_2_BC)
            
            # STEP 5.2: MITOTIC Cell Division for DIPLOID CELLS
            
            # STEP 5.21: Prepare the data types and format to send it to CELL DIVISION MODULE            
            # Grouping the cells based on their cell division time for "DIFFERENTIAL REPRODUCTION".
            cell_groups =  Yeast_Simulator_Subprograms.group_cells(evolved_cells,n_evolved_cells)

            # Decide the parameters
            n_BC_generations = 5
            N_BN = 1           
                        
            # Expanding the evolved_cells size
            genotypes_1 = BackCrossing.expand_format(evolved_cells,n_BC_generations,N_BN,n_evolved_cells)
            
            # RE-Opening the HAPLOTYPE files to append the new HAPLOTYPE OF POINT MUTATIONS AND GENE DUPLICATIONS.
            point_mutation = open("Point_Mutations.txt","a")
            
            # Assymetrical Cell division for "DIPLOIDS"
            ploidy = 2; div_type = 2; s_size = n_evolved_cells; fix_cal = 0            
            (genotypes_1,PM_2_BC,PM_fitness_2_BC,PM_strand,n_mut) = Yeast_Simulator_Cell_Division.asymmetrical_cell_division(genotypes_1,cell_groups,PM_fitness_2_BC,PM_temp,PM_strand,n_BC_generations,N_BN,n_evolved_cells,n_mut,point_mutation,ploidy,"DBC",div_type,s_size,fix_cal,PMT)
            
            #for ii in genotypes_1:
            #    print "After:",ii
            
            # Closing the HAPLOTYPE files
            point_mutation.close()
            
            # Opening the Point mutation Haplotype files
            f1 = open("Point_Mutations.txt","r")
            PM_hap_2_BC = f1.readlines()
            f1.close()
            
            #LOST3 = BackCrossing.freq_test3(PM_2_BC,initial_PMs,550,PM_fitness_2_BC)
            #print "LOST3:",i,"\t",j,"\t",LOST3            
            
            #######################################################################################################
            # STEP VI: Sub Sample 50000 of them: 
            #######################################################################################################            
            
            # Selecting the highly evolved cells from the population to backcross with founder population.
            # But here we use this funtion just to convert the format
            # evolved_cells,GD_2_BC,GD_fitness_2_BC,GD_hap_2_BC,PM_2_BC,PM_fitness_2_BC,PM_hap_2_BC = BackCrossing.select_evolved_cells_randomly_to_cross(n_evolved_cells,genotypes_1,n_cell,GD_2_BC,GD_fitness_2_BC,GD_hap_2_BC,PM_2_BC,PM_fitness_2_BC,PM_hap_2_BC)
            
            # GD_2_BC = BackCrossing.convert_format_1(GD_2_BC)
            PM_2_BC = BackCrossing.convert_format_1(PM_2_BC)
            
            #######################################################################################################
            # STEP VII: Force 99% of diploids into haploids
            #######################################################################################################
            
            (genotypes,PM,PM_fitness) = BackCrossing.diploids_2_haploids(genotypes_1,PM_2_BC,PM_fitness_2_BC,PM_strand,s_size)
            
            #(genotypes,GD,GD_fitness,PM,PM_fitness) = BackCrossing.diploids_2_haploids_1(genotypes_1,GD_2_BC,GD_fitness_2_BC,GD_strand,PM_2_BC,PM_fitness_2_BC,PM_strand,n_evolved_cells)
            #######################################################################################################
            # STEP VIII: Sub Sample 50000 individuals from population: THIS STEP CAN BE AVOIDED - code rewiring needed
            #######################################################################################################
            (genotypes,PM,n_evolved_cells) = BackCrossing.sampling_cells(genotypes,PM,PM_fitness,n_evolved_cells)
            # FREQ3 = BackCrossing.freq_test1(PM,500)
            # for ii in genotypes:
            #    if ii[5] ==2:
            #        print "***:",ii
            #    else:
            #        print ii
            #sys.exit()
            #######################################################################################################
            # STEP IX: Let Haploids reproduce mitotically for 5 Gs not diploids
            # STEP X: Sub sample 50000 from population
            # STEP XI: OPTIONAL: Repeat step 9 and 10
            #######################################################################################################
            # Mitotic division Parameters for Non-selection exp
            ploidy = 1; div_type = 1; fix_cal = 0            
            # No. of non-selection rounds
            NR_NS = 1
            # No. of selection rounds
            NR_SE = 5
            total_rounds = NR_NS + NR_SE
            
            for k in range(total_rounds):
                ###################################################################################################
                # STEP XII, STEP XIII and STEP XIV: Allow mitotic cell division "UNDER STRESS" for N generations
                # Instead of writing extra code for SELECTION, selection part is included in this LOOP.
                # IF CONDITION will decide when the SELCTION EXPERIMENT will be started after NON-SELCTION EXPs
                ###################################################################################################
                
                if k == NR_NS:
                    genotypes = BackCrossing.prepare_for_s_exp(genotypes,PM,PM_fitness,PM_strand,n_evolved_cells)

                    # Mitotic division Parameters for selection exp
                    ploidy = 1; div_type = 1; fix_cal = 0                    
                
                # Grouping the cells based on their cell division time for "DIFFERENTIAL REPRODUCTION".
                cell_groups =  Yeast_Simulator_Subprograms.group_cells(genotypes,n_evolved_cells)

                # Decide the parameters
                
                #if k == NR_NS:
                #    n_BC_generations = 5
                #    N_BN = 5           
                #else:
                #    n_BC_generations = 5
                #    N_BN = 1
                
                n_BC_generations = 5
                N_BN = 1                    
                
                # RE-Opening the HAPLOTYPE files to append the new HAPLOTYPE OF POINT MUTATIONS AND GENE DUPLICATIONS.
                point_mutation = open("Point_Mutations.txt","a")
                
                # Higher sample size is needed for backcrossing module
                if k == total_rounds-1:
                    s_size = n_evolved_cells * 2
                else:
                    s_size = n_evolved_cells

                # Mitotic Cell Division: Assymetrical Cell division for "Haploids"
                (genotypes,PM,PM_fitness,PM_strand,n_mut) =  Yeast_Simulator_Cell_Division.asymmetrical_cell_division(genotypes,cell_groups,PM_fitness,PM,PM_strand,n_BC_generations,N_BN,n_evolved_cells,n_mut,point_mutation,ploidy,"DBC",div_type,s_size,fix_cal,PMT)
                
                # This one calculates the Haplotype frequencies after BC
                prefered_BCs = [4,9]
                if k == total_rounds-1 and i in prefered_BCs:
                    Yeast_Simulator_Cell_Division.calculate_haplotype_freq_1(genotypes,PM,n_evolved_cells,i)
                    
                # Closing the HAPLOTYPE files
                point_mutation.close()
                # gene_dup.close()
            print "END_OF_SELECTION_EXP"
            #LOST5 = BackCrossing.freq_test3(PM,initial_PMs,500,PM_fitness_2_BC)
            #print "LOST5:",i,"\t",j,"\t",LOST5
            # END OF STEP 9-10 REPEAT            
            #######################################################################################################
            # STEP XII: Allow mitotic cell division under stress for 5 generations
            #######################################################################################################
            
            #######################################################################################################
            # STEP XIII: Sub Sample 50000 cells
            #######################################################################################################
            
            #######################################################################################################
            # STEP XIV: OPTIONAL: Repeat 12-13 
            #######################################################################################################
            
            #######################################################################################################
            # STEP XV: Take sub sample to STEP 1
            #######################################################################################################
            
            #######################################################################################################
            # STEP XVI: Repeat 1-16 steps number of times
            #######################################################################################################

            # Find the Fitnesss type
            # fitness_type = BackCrossing.find_fitness_type(GD,GD_fitness)
            
            # Calculate allele frequency in the selected cells.
            # BackCrossing.testing_baba(GD_2_BC,GD_fitness_2_BC,4)
            
            #GDs,GD_new_ids = BackCrossing.combine_fitness_1(GDs,GD,GD_fitness)
            #PMs,PM_new_ids = BackCrossing.combine_fitness_1(PMs,PM,PM_fitness)
                   
            # initial_GDs_temp = GD_fitness.keys()
            initial_PMs_temp = PM_fitness.keys()
            text = "ABC" + str(i+1)
            # Preferred BC results
            ABCs = [1,3,5,8,10]
            if i+1 in ABCs:
                #print text,i
                #BackCrossing.calculate_frequency('GD',text,j,i+1,GD,initial_GDs,GD_fitness,param_text)
                BackCrossing.calculate_frequency('PM',text,j,i+1,PM,initial_PMs,PM_fitness,param_text,s_size)
            print "###################################################"


            #GD_2_BC = BackCrossing.convert_format_1(GD)
            PM_2_BC = BackCrossing.convert_format_1(PM)
            
            # Opening the Point mutation Haplotype files
            f1 = open("Point_Mutations.txt","r")
            PM_hap_2_BC = f1.readlines()
            f1.close()
                            
            PM_hap_2_BC = BackCrossing.convert_format_2(PM_2_BC,PM_hap_2_BC)
            #GD_hap_2_BC = BackCrossing.convert_format_2(GD_2_BC,GD_hap_2_BC)
            
            evolved_cells = {}; n_dip = 0
            for kk in range(n_evolved_cells*2):
                if genotypes[kk][5] == 2:
                    n_dip += 1
                #evolved_cells[kk] = deepcopy(genotypes[kk])
                evolved_cells[kk] = genotypes[kk]
                
            print "Diploids:",n_dip
            print "END_OF_ABC"
        #print "TIME AFTER:",j,datetime.datetime.now() 
        print "END_OF_BC_REPEAT"
    #epistasis.close()
    PMT.close()
#############################################################################################################
#############################################################################################################

#cProfile.run("Yeast_lab()")


# This is to run this module normally like "python Yeast_Simulator.py"
if __name__ == "__main__":
    import sys
    Yeast_lab()
