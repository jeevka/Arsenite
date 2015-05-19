from __future__ import division
import sys
import random
import math
import numpy
import scipy.stats

############################################################################################################
# Importing user defined modules.
############################################################################################################
import Yeast_Simulator
import Yeast_Simulator_Mutations as PM
import Yeast_Simulator_Subprograms as SP
import Artificial_Cells
############################################################################################################
# Global values - constants
############################################################################################################
gene_duplication_fitness = Yeast_Simulator.gene_duplication_fitness
shape = Yeast_Simulator.gene_duplication_shape
scale = Yeast_Simulator.gene_duplication_scale
#shape = 2
#scale = 1
#gene_duplication_fitness = Yeast_Simulator.gene_duplication_random_fitness
#GD_haplotype = Yeast_Simulator.gene_duplication_haplotype
beneficial_rate = float(Yeast_Simulator.gene_duplication_beneficial/100)
neutral_rate = float(Yeast_Simulator.gene_duplication_neutral/100)
deleterious_rate = float(Yeast_Simulator.gene_duplication_deleterious/100)
fitness_truncation = Yeast_Simulator.gene_duplication_fitness_truncation
############################################################################################################

############################################################################################################
# This function will introduce the gene duplication either to mother or daughter cell
############################################################################################################
def gene_duplications(mother_cell,daughter_cell,k,tp,n_gene_dup,gene_dup,GD,GD_fitness,GD_strand,ploidy,fix_cal):
    
    # Fitness
    # fix_cal decides whether GDs carry any fitness effects.
    if fix_cal == 0:
        fitness = 0
    else:    
        fitness = gene_duplication_fitness[n_gene_dup]
    
    # Haplotype
    # haplotype = GD_haplotype[n_gene_dup]
            
    # Choosing between Mother and daughter cells.    
    # Mother cell
    if random.random() >= 0.5:
        # Calcualting the new cell division time based on new fitness.
        if fitness != 0:
            mother_cell[1] = SP.calculate_cell_division_time(mother_cell[1],fitness)
        # Assigning the Geneduplication Structure.
        gene_dup.write(str(PM.assign_haplotype(mother_cell[5])))
        gene_dup.write("\t")

        # Adding the gene duplication fitness
        GD_fitness[n_gene_dup] = fitness
        GD_strand[n_gene_dup] = random.randrange(1,ploidy+1)
        mother_cell[4] = mother_cell[4] + 1
        GD[k].append(n_gene_dup)

    # Daughter cell
    else:
        # Calcualting the new cell division time based on new fitness.
        if fitness != 0:
            daughter_cell[1] = SP.calculate_cell_division_time(daughter_cell[1],fitness)

        #haplotype = decide_chr_type(daughter_cell[5],haplotype)    
        gene_dup.write(str(PM.assign_haplotype(daughter_cell[5])))
        gene_dup.write("\t")
        
        # Adding the gene duplication fitness
        #GD_fitness.append(fitness)
        GD_fitness[n_gene_dup] = fitness
        GD_strand[n_gene_dup] = random.randrange(1,ploidy+1)

        daughter_cell[4] = daughter_cell[4] + 1
        GD[tp].append(n_gene_dup)
    gene_dup.write("\n")

    return mother_cell,daughter_cell,n_gene_dup,GD,GD_fitness,GD_strand

##############################################################################################################
# Decide whether you choose the homolog chr or not
##############################################################################################################
def decide_chr_type(type,haplotype):
    print type
    print haplotype
    sys.exit()

##############################################################################################################
# Calulate the fitness level
##############################################################################################################
def calculate_fitness(truncation):
    ran_n = random.random()
    fitness_limit = float(50)/float(truncation)
    fitness = float(ran_n)/float(fitness_limit)
    
    return fitness
