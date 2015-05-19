from __future__ import division
import sys
import random
import operator
from copy import copy
from copy import deepcopy
import math
import re
import time
from datetime import timedelta
import datetime
#from numpy import random, array, vstack, column_stack, ones, zeros
from Numeric import *
 
import profile
import pstats

from decimal import *
import meiosis
import Yeast_Model_SubPrograms_Cython



# Global variables
haploid_cells = {}
n_hap = 0
chr_phy_map = {}
backcross_haploid_cells = {}
minimum_cell_division_time = 0
maximum_cell_division_time = 0

###################################################################################################
# Here we choose the ten best cells with faster cell division time.
###################################################################################################
def select_evolved_cells_to_cross(number_of_cells, genotypes, minimum_cell_division_time,maximum_cell_division_time):
    #minimum_cell_division_time = 90
    #maximum_cell_division_time = 900
    selected_cells = {}
    k = 0
    for i in range(1,len(genotypes)+1):
        if genotypes[i][2] <= maximum_cell_division_time:
           maximum_cell_division_time = genotypes[i][2]
           number_of_cells =  number_of_cells - 1
           selected_cells[k] = genotypes[i]
           k = k + 1
           if number_of_cells == 0:
              break
                
    return selected_cells        

###################################################################################################
# Here we allow the wild and evolved cells to mate
###################################################################################################
def back_cross(wild_type_individuals,evolved_cells,mutations,beneficial_mutations,beneficial_mutations_fitness,gene_deletions,gene_deletions_fitness,gene_duplications,gene_duplications_fitness):
    diploid_cells = {}
    # Centimorgan values for yeast chrs. Ref: www.yeastgenome.org/pgMaps/pgMap.shtml
    cMorgan_values = [0.45,0.30,0.48,0.31,0.35,0.48,0.37,0.30,0.45,0.30,0.38,0.36,0.34,0.37,0.33,0.29]
    number_of_cells = 0
    for i in range(0,len(evolved_cells)):
        number_of_cells = number_of_cells + 4
        # Build chr physical map by combining wild and evolved ones.
        diploid_cells[i] = build_chromosome_physical_map(wild_type_individuals[i],evolved_cells[i],gene_deletions,gene_duplications,mutations)
        # print diploid_cells[i]
        # calculate rebination rate
        # diploid_cells[i] = calculate_recombine_rate(diploid_cells[i],cMorgan_values)
        # Once the recombination is done. Lets allow the DNA to replicate.
        dna_replication_and_division_one(diploid_cells[i],cMorgan_values)
        # After recombination each and every mutation will be shuffled between chromosomes.
        # So based on the beneficial mutation we have to calculate the cell division time.
        cell_division_time = calculate_cell_division_time(number_of_cells,wild_type_individuals[i],evolved_cells[i],mutations,beneficial_mutations,beneficial_mutations_fitness,gene_deletions,gene_deletions_fitness,gene_duplications,gene_duplications_fitness)
        # Create the cells Based on the Genetic variations and cell division time.
        # create_cells(cell_division_time)
    
    return backcross_haploid_cells


###################################################################################################        
# This will build the physical map of chromosomes with mutations
###################################################################################################
def build_chromosome_physical_map(wild_type,evolved_cell,gene_deletions,gene_duplications,mutations):
    # This array of array(list of lists) contains 32+1 arrays inside.32 for diploid number of Yeast chromosomes. And the
    # Last array is to store the cell division time.
    chr_phy_map[0] = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    id = ""
    # building map by Mutations
    for i in range(0,len(evolved_cell[7])):
        # its so foolish to do in this way. Have to change it Soon.
        for j in range(0,len(mutations)):
            if evolved_cell[7][i] == mutations[j]:
               id = j         
        (chr_no, chr_position) = split_chro_position(evolved_cell[7][i])
        t = str(0.1) + str(id)
        chr_position = float(chr_position) +  float(t)
        s = len(str(id))
        d = "%."+ str(s) + "f"
        chr_position = d % chr_position
        chr_phy_map[0][chr_no*2 -1].append(float(chr_position))
        chr_phy_map[0][chr_no*2 -2].append(0)
        chr_phy_map[0][chr_no*2 -1].sort()
 
    # Building map by gene duplications
    for i in range(0,len(evolved_cell[8])):
        id = evolved_cell[8][i]
        (chr_no, chr_position) = split_chro_position(gene_duplications[id])
        t = str(0.2) + str(id)
        chr_position = float(chr_position) +  float(t)
        s = len(str(id))
        d = "%."+ str(s) + "f"
        chr_position = d % chr_position
        chr_phy_map[0][chr_no*2 -1].append(float(chr_position))
        chr_phy_map[0][chr_no*2 -2].append(0)
        chr_phy_map[0][chr_no*2 -1].sort()
    
    # Building map by gene deletions
    for i in range(0,len(evolved_cell[10])):
        id = evolved_cell[10][i]
        (chr_no, chr_position) = split_chro_position(gene_deletions[id])
        t = str(0.3) + str(id)
        chr_position = float(chr_position) +  float(t)
        s = len(str(id))
        d = "%."+ str(s) + "f"
        chr_position = d % chr_position
        chr_phy_map[0][chr_no*2 -1].append(float(chr_position))
        kl = len(chr_phy_map[0][chr_no*2 -1])
        chr_phy_map[0][chr_no*2 -2].append(0)
        chr_phy_map[0][chr_no*2 -1].sort()

    return chr_phy_map

###################################################################################################
# This will calculate the physical distance between markers in the chr
# and using that it will calculate recombination fractions and do 
# recombination between homologous chromosomes.
###################################################################################################
def calculate_recombine_rate(diploid_cell,cMorgan_values):
    physical_distance_list = []        
    for j in range(0,16):
        # Will tell you whether there is a mutation or not.
        mutation_indicator = 0
        if len(diploid_cell[0][j*2+1]) > 1:
            mutation_indicator = 1
            # if the chr has zero or one mutation
            if diploid_cell[0][j*2+1].count(0) >= len(diploid_cell[0][j*2+1])-1:
                continue
            # if the chr has more than one mutation
            else:
                diploid_cell[0][j*2+1].sort()
                for k in range(0,len(diploid_cell[0][j*2+1])-1):
                   # finding the physical distance and dividing by 1000 to convert it to cM/Kpbs
                   # and multiplying by the cM value
                   if diploid_cell[0][j*2+1][k+1] != 0 and diploid_cell[0][j*2+1][k] != 0:
                    distance = (((diploid_cell[0][j*2+1][k+1]) - (diploid_cell[0][j*2+1][k]))/1000) * cMorgan_values[j]
                    physical_distance_list.append(distance)
                    # Calling invHaldane to convert map distances recombination fractions
                    n = len(physical_distance_list)
                    physical_distance_list[n-1] = meiosis.invHaldane(physical_distance_list[n-1])

        # Converting the chr map to two dimensional array
        # So that we can send this to arne's meiosis program to recombine.
        # if there is more than one marker
        if mutation_indicator == 1:
            n_markers = len(diploid_cell[0][j*2+1])
            chr = array([0.0,0.0])
            diploid_chr_map = array([resize(chr,(n_markers,2))])
            for k in range(0,len(diploid_cell[0][j*2+1])):
                diploid_chr_map[0][k][1] = diploid_cell[0][j*2+1][k]
                diploid_chr_map[0][k][0] = diploid_cell[0][j*2][k]

            diploid_chr_map[0] = meiosis.recombine(diploid_chr_map[0],physical_distance_list)
            # repack them
            for k in range(0,len(diploid_cell[0][j*2+1])):
                diploid_cell[0][j*2+1][k] = diploid_chr_map[0][k][1] 
                diploid_cell[0][j*2][k] = diploid_chr_map[0][k][0] 
                
    return diploid_cell

###################################################################################################
# Here the DNA replicates and the cells divide twice to make 4 haploids.
###################################################################################################
def dna_replication_and_division_one(diploid_cell,cMorgan_values):
    diploid_1 = []
    diploid_2 = []
    # DNA replicate
    from copy import copy
    chr_copy1 = copy(diploid_cell)
    chr_copy2 = copy(diploid_cell)
    # It will calculate the recombination rate based on the genetic distance and
    # cMorgan values. And then it will call recombine funtion and return the
    # recomnbined chrmosomes. 
    chr_copy1 = calculate_recombine_rate(chr_copy1,cMorgan_values)
    chr_copy2 = calculate_recombine_rate(chr_copy2,cMorgan_values)
    
    # Forming two diploid cells
    # Diploid cell one.
    for i in range(0,16): 
        # choosing the left chr for each pair
        if random.randrange(0,2) == 0:
            diploid_1.append(chr_copy1[0][i*2])
            diploid_2.append(chr_copy2[0][i*2])
        else:
            diploid_1.append(chr_copy2[0][i*2])
            diploid_2.append(chr_copy1[0][i*2])
        # choosing the right chr for each pair
        if random.randrange(0,2) == 0:
            diploid_1.append(chr_copy1[0][i*2+1])
            diploid_2.append(chr_copy2[0][i*2+1])
        else:
            diploid_1.append(chr_copy2[0][i*2+1])
            diploid_2.append(chr_copy1[0][i*2+1])
    
    # Calling the function for cell division to become 4 haploid cells.
    diploids_to_haploids(diploid_1,diploid_2)
    # Here we consider replicate_dna_copy as a one diploid cell and
    # diploid_cell as another diploid cell from the parent.
    # cell division-I
    # for i in range(0,16):
        # Each chromosome is replicated. So for each chrosome
        # We should choose the copy randomly. It will produce two cells
        # Including the replication copy, there are four copies of each chr strand.
        # We have to group them as two pairs randomly.
        # 1-wildtype original, 2-evolved original, 3-wiltype copy,4-evolved copy

###################################################################################################
# Here we seperate the two diploid cell divides into four haploid cells.
###################################################################################################
def diploids_to_haploids(diploid_1,diploid_2):
    hap1 = [];hap2 = [];hap3 =[]; hap4 = []
    n_hap = len(haploid_cells)
    for i in range(0,16):
        if random.randrange(0,2) == 0: 
           hap1.append(diploid_1[i*2])
           hap2.append(diploid_1[i*2+1])
           hap3.append(diploid_2[i*2])
           hap4.append(diploid_2[i*2+1])
        else:
           hap1.append(diploid_1[i*2+1])
           hap2.append(diploid_1[i*2])
           hap3.append(diploid_2[i*2+1])
           hap4.append(diploid_2[i*2])          
    haploid_cells[n_hap+1] = hap1
    haploid_cells[n_hap+2] = hap2
    haploid_cells[n_hap+3] = hap3
    haploid_cells[n_hap+4] = hap4


###################################################################################################
# Here we calculate the cell division time and build the cell based on all the Genetic variations.
###################################################################################################
def calculate_cell_division_time(number_of_cells,wild_type_individuals,evolved_cells,mutations,beneficial_mutations,beneficial_mutations_fitness,gene_deletions,gene_deletions_fitness,gene_duplications,gene_duplications_fitness):
    n_hap = len(haploid_cells)
    cell_division_time = []
    minimum_cell_division_time = 90
    maximum_cell_division_time = 900
    m = 0
    s = number_of_cells-3
    e = number_of_cells+1
    for i in range(s,e):
        m = m + 1
        # Building the cells
        backcross_haploid_cells[i] = [0,900,900,16,0,0,0,{},[],[],[],0]
        cell_division_time.append(maximum_cell_division_time)
        for j in range(0,16):
            for k in range(0,len(haploid_cells[i][j])):
                temp_1 = re.split("\.",str(haploid_cells[i][j][k]))
                chr_pos = int(temp_1[0])
                if len(temp_1) > 1:
                    temp_2 = str(temp_1[1])
                    genetic_variation_id = temp_2[0]
                    
                    # Mutations and fitness
                    if int(genetic_variation_id) == 1:
                        id = ""
                        for l in range(1,len(temp_2)):
                            id = id + temp_2[l]
                        # adding the mutations to the cells.
                        l = len(backcross_haploid_cells[i][7])
                        if l == 0:
                            l = 0
                        else:
                            l = l -1
                        if id == "":
                            id = 0

                        # Adding Mutations into the Cell.    
                        backcross_haploid_cells[i][7][l] = mutations[int(id)] 
                        # Get the mutation ID.
                        for l in range(0,len(beneficial_mutations)):
                            if mutations[int(id)] == beneficial_mutations[l]:
                                cell_division_time[m-1] = cell_division_time[m-1] - ((cell_division_time[m-1] - minimum_cell_division_time) * beneficial_mutations_fitness[l])
                                # Adding the beneficial mutation id to the cell.
                                backcross_haploid_cells[i][10].append(l)
                        #print "2:",cell_division_time[m-1]
                        # Re-calculating the cell divisioin time based on beneficial mutations.
                        backcross_haploid_cells[i][2] = cell_division_time[m-1]
                    
                    # Gene duplications and fitness.
                    if int(genetic_variation_id) == 2:
                        id = ""
                        for l in range(1,len(temp_2)):
                            id = id + temp_2[l]
                        if id == "":
                            id = 0
                        # Adding the Gene duplications to the Cell.    
                        backcross_haploid_cells[i][8].append(int(id))     
                        # Calcualting the cell division time from Gene duplications.
                        cell_division_time[m-1] = cell_division_time[m-1] - ((cell_division_time[m-1] - minimum_cell_division_time) * gene_duplications_fitness[int(id)])    
                        #print "4:",cell_division_time[m-1]
                        # Assigning the cell new cell division time to the cell.
                        backcross_haploid_cells[i][2] = cell_division_time[m-1]
                        
                    
                    # Gene deletions and fitness    
                    if int(genetic_variation_id) == 3:
                        id = ""
                        for l in range(1,len(temp_2)):
                            id = id + temp_2[l]
                        if id == "":
                            id = 0
                        # Adding the Gene deletions to the cell.    
                        backcross_haploid_cells[i][9].append(int(id))
                        # Calcualting the cell division time from Gene deletion.
                        cell_division_time[m-1] = cell_division_time[m-1] - gene_deletions_fitness[int(id)]
                        
                        # Assigning new cell division time.
                        backcross_haploid_cells[i][2] = math.ceil(cell_division_time[m-1])
        # Assigning the next cell division time
        backcross_haploid_cells[i][1] = math.ceil(backcross_haploid_cells[i][2])
    return cell_division_time    

###################################################################################################
# Each dict contains all the information for cells.
# We need to split the chr no and position from that.
###################################################################################################
def split_chro_position(haplotype):
    temp_list = re.split(":",str(haplotype))
    temp_var = re.split("{",temp_list[0])
    chr_no = int(temp_var[1])
    chr_position = re.split("}",temp_list[1])    
    
    return chr_no,chr_position[0]
