from __future__ import division
from datetime import timedelta
from copy import copy
from copy import deepcopy
from numpy import *
from collections import Counter
import sys
import random
import operator
import math
import re
import time
import datetime
import profile
import pstats
import textwrap

################################# User defind modules ######################################################
#import Yeast_Simulator
import meiosis
#import Yeast_Simulator_Subprograms as YSP
#import Yeast_Simulator_Cell_Division

##########################################################################################################
# Global variables
##########################################################################################################

##########################################################################################################
# Choose the N number of wildtype cells before the experiment to backcross
##########################################################################################################
def create_wild_type_cells():
    WT_cells = [0,0,0,0,0,0] #ones((1,6),dtype=int32)
    
    WT_cells[0] = 180
    WT_cells[1] = 180
    WT_cells[2] = 16
    WT_cells[5] = 1
	
    return WT_cells

##########################################################################################################
# Here we choose the number of best cells with faster cell division time.
##########################################################################################################
def select_evolved_cells_to_cross_1(number_of_cells, genotypes,n_individuals,PM,PM_fitness,PM_hap):    
    
    # To choose the best cells or highly adapted cells.
    best_CD_time = find_best_cells(genotypes,n_individuals)
    
    selected_cells = ones((number_of_cells,6),dtype=int32)
    k = 0
    
    s_PM = []; s_PM_fitness = []; s_PM_hap = []
    
    # Initializing the lists
    for i in xrange(number_of_cells):
	s_PM_fitness.append([]);s_PM_hap.append([])
    
    # Choosing the cells and their Haplotypes with fitness
    for j in xrange(len(best_CD_time)):
	for i in xrange(1,n_individuals):
	    if genotypes[i][1] <= best_CD_time[j]:
	       time = genotypes[i][1]
	       number_of_cells =  number_of_cells - 1
	       selected_cells[k] = genotypes[i]
	       selected_cells[k][0] = genotypes[i][1]

	       # Copying the cell's GD ids
	       #s_GD.append(GD[i])	       
	       #s_GD_fitness,s_GD_hap = copy_cat(GD[i],GD_fitness,GD_hap,s_GD,s_GD_fitness,s_GD_hap,k)
	       
	       # Copying cell's PM ids
	       s_PM.append(PM[i])
	       s_PM_fitness,s_PM_hap = copy_cat(PM[i],PM_fitness,PM_hap,s_PM,s_PM_fitness,s_PM_hap,k)
	       
	       k += 1
	       if number_of_cells == 0:
		  break
		
	if number_of_cells == 0:
	    break
	
    return selected_cells,s_PM,s_PM_fitness,s_PM_hap        


##########################################################################################################
# Here we choose the number of best cells with faster cell division time.
##########################################################################################################
def select_evolved_cells_to_cross(number_of_cells, genotypes,n_individuals,PM,PM_fitness,PM_hap):    
    
    # To choose the best cells or highly adapted cells.
    best_CD_time = find_best_cells(genotypes,n_individuals)
    
    selected_cells = ones((number_of_cells,6),dtype=int32)
    k = 0
    
    s_PM = []; s_PM_fitness = {}; s_PM_hap = {}
    
    # Choosing the cells and their Haplotypes with fitness
    for j in xrange(len(best_CD_time)):
	for i in xrange(1,n_individuals):
	    if genotypes[i][1] <= best_CD_time[j]:
	       time = genotypes[i][1]
	       number_of_cells =  number_of_cells - 1
	       selected_cells[k] = genotypes[i]
	       selected_cells[k][0] = genotypes[i][1]

	       # Copying the cell's GD ids
	       #s_GD.append(GD[i])	       
	       #s_GD_fitness,s_GD_hap = copy_cat(GD[i],GD_fitness,GD_hap,s_GD,s_GD_fitness,s_GD_hap,k)
	       
	       # Copying cell's PM ids
	       s_PM.append(PM[i])
	       s_PM_fitness,s_PM_hap = copy_cat(PM[i],PM_fitness,PM_hap,s_PM,s_PM_fitness,s_PM_hap,k)
	       
	       k += 1
	       if number_of_cells == 0:
		  break
		
	if number_of_cells == 0:
	    break
	
    return selected_cells,s_PM,s_PM_fitness,s_PM_hap        

##########################################################################################################
# Here we choose the number of cells randomly
##########################################################################################################
def select_evolved_cells_randomly_to_cross(number_of_cells, genotypes,n_individuals,PM,PM_fitness,PM_hap):
    
    choosen_cells = random.sample(xrange(n_individuals),number_of_cells)
 			     
    selected_cells = ones((number_of_cells,6),dtype=int32)

    #s_GD = []; s_GD_fitness = {}; s_GD_hap = {}
    s_PM = []; s_PM_fitness = {}; s_PM_hap = {}

    
    # Choosing the cells and their Haplotypes with fitness
    for i in xrange(number_of_cells):
	time = genotypes[i][1]
	id = choosen_cells[i]

	selected_cells[i] = genotypes[id]
	selected_cells[i][0] = genotypes[id][1]

	# Copying the cell's GD ids
	#s_GD.append(GD[id])
	#s_GD_fitness,s_GD_hap = copy_cat(GD[id],GD_fitness,GD_hap,s_GD,s_GD_fitness,s_GD_hap,i)
	       
	# Copying cell's PM ids
	s_PM.append(PM[id])
	s_PM_fitness,s_PM_hap = copy_cat(PM[id],PM_fitness,PM_hap,s_PM,s_PM_fitness,s_PM_hap,i)
	
    return selected_cells,s_PM,s_PM_fitness,s_PM_hap

##########################################################################################################
# Here we choose the number of cells randomly
##########################################################################################################
def select_evolved_cells_randomly_to_cross_1(number_of_cells, genotypes,n_individuals,GD,GD_fitness,GD_hap,PM,PM_fitness,PM_hap):
    
    choosen_cells = random.sample(xrange(n_individuals),number_of_cells)
 			     
    selected_cells = ones((number_of_cells,6),dtype=int32)

    s_GD = []; s_GD_fitness = []; s_GD_hap = []
    s_PM = []; s_PM_fitness = []; s_PM_hap = []

    # Initializing the lists
    for i in xrange(number_of_cells):
        s_GD_fitness.append([]);s_GD_hap.append([])
	s_PM_fitness.append([]);s_PM_hap.append([])
    
    # Choosing the cells and their Haplotypes with fitness
    for i in xrange(number_of_cells):
	time = genotypes[i][1]
	id = choosen_cells[i]

	selected_cells[i] = genotypes[id]
	selected_cells[i][0] = genotypes[id][1]

	# Copying the cell's GD ids
	s_GD.append(GD[id])
	s_GD_fitness,s_GD_hap = copy_cat(GD[id],GD_fitness,GD_hap,s_GD,s_GD_fitness,s_GD_hap,i)
	       
	# Copying cell's PM ids
	s_PM.append(PM[id])
	s_PM_fitness,s_PM_hap = copy_cat(PM[id],PM_fitness,PM_hap,s_PM,s_PM_fitness,s_PM_hap,i)

    return selected_cells,s_GD,s_GD_fitness,s_GD_hap,s_PM,s_PM_fitness,s_PM_hap
##########################################################################################################
# Select only the DIPLOID CELLS
##########################################################################################################
def select_diploid_cells(genotypes,PM,PM_fitness,PM_hap):

    # To choose the best cells or highly adapted cells.
    import random
    choosen_cells = []
    number_of_cells = 0
    for i in xrange(len(genotypes)):
	if genotypes[i][5] == 2:
	    number_of_cells += 1
	    choosen_cells.append(i)
	    
    selected_cells = ones((number_of_cells,6),dtype=int32)
    
    # s_GD = []; s_GD_fitness = []; s_GD_hap = []
    s_PM = []; s_PM_fitness = []; s_PM_hap = []
    
    # Choosing only the Diploid cells
    j = 0
    for i in xrange(len(choosen_cells)):
	id = choosen_cells[i]
    
	selected_cells[i] = genotypes[id]
	selected_cells[i][0] = genotypes[id][1]
    
	# Copying the cell's GD ids
	#s_GD.append(GD[id])
	#s_GD_fitness.append(GD_fitness[id])
	#s_GD_hap.append(GD_hap[id])
	
		   
	# Copying cell's PM ids
	s_PM.append(PM[id])
	s_PM_fitness.append(PM_fitness[id]);
	s_PM_hap.append(PM_hap[id])
	
    return selected_cells,s_PM,s_PM_fitness,s_PM_hap,number_of_cells

##########################################################################################################
# Just to copy
##########################################################################################################
def copy_cat(ID,fitness,hap,s_ID,s_fitness,s_hap,L1):
    for i in ID:
	if s_fitness.has_key(i):
	    pass
	else:
	    s_fitness[i] = fitness[i]
	    s_hap[i] = hap[i].strip()
	    
	"""
	s_fitness[L1].append(fitness[i])
	s_hap[L1].append(hap[i].strip())
	"""
	
    return s_fitness,s_hap

##########################################################################################################
# Just to copy
##########################################################################################################
def copy_cat_1(ID,fitness,hap,s_ID,s_fitness,s_hap,L1):
    for i in xrange(len(ID)):
        id = ID[i]
	s_fitness[L1].append(fitness[id])
	s_hap[L1].append(hap[id].strip())
	
    return s_fitness,s_hap

##########################################################################################################
# Sort the cell's cell division time to find the top best cells based on their cell division times
##########################################################################################################
def find_best_cells(genotypes, n_individuals):
    CDT = {}
    for i in xrange(n_individuals):
	time = genotypes[i][1]
	CDT[time] = time
	
    return CDT.keys()
	
##########################################################################################################	
# Here we allow the wild and evolved cells to mate
##########################################################################################################
def back_cross(evolved_cells,PM,PM_fitness,PM_hap,initial_PMs,j,i,p_m,PM_strand):
    
    # Calculate Mating cells proportion from p_m
    n1,n2,n3,n4 = calculate_mating_params(evolved_cells,p_m)
    
    selected_cells = random.sample(range(0,len(evolved_cells)),len(evolved_cells))
    n = len(evolved_cells) * 2**5
    s_cells = []
    genotypes = ones((n,6),dtype=int32)
    
    diploid_cells = {}
        
    # To store the Gene duplications and Mutations
    PM_new = []; 
    
    # Centimorgan values for yeast chrs. Ref: www.yeastgenome.org/pgMaps/pgMap.shtml
    cMorgan_values = [0.45,0.30,0.48,0.31,0.35,0.48,0.37,0.30,0.45,0.30,0.38,0.36,0.34,0.37,0.33,0.29]
    n_cell = 0
    k = 0
    
    # only 50% of the cells mates
    n_evolved_cells = int(n4) #int(len(evolved_cells)/2)
    
    # Number of evolved cells
    n_e_c = 0;
    # Number of Diploid cells 
    n_d_c = 0
    wt = 0
    nn1 = nn2 = nn3 = 0
    for i in xrange(n_evolved_cells):
	# Randomly decide whether the cells are going to mate and form Diploids or not !
	mate_indicator = random.random()
	
	# Evolved Haploids
	if mate_indicator <= n3:
	    sc = selected_cells[n_e_c]
	    genotypes,PM_new = copy_sub(n_cell,genotypes,evolved_cells[sc],PM_new,PM[sc])
	    n_e_c += 1
	    n_cell += 1
	    
	    # Storing the selected_cells to calculath PM, GD freq
	    s_cells.append(k)
	    
	    k += 1
	    nn1 += 1
	# Wild Type Haploids
	elif mate_indicator >n3 and mate_indicator <= n2:
	    wt += 1
	    genotypes[n_cell] = create_wild_type_cells()
	    PM_new.append([])
	    
	    n_cell += 1
	    k += 1
	    nn2 += 1
	# Mating between Haploids and formation of Diploids	    
	else:
	    nn3 += 1
	    n_e_c += 1
	    sc = selected_cells[n_e_c]
	    
	    # Storing the selected_cells to calculath PM, GD freq
	    s_cells.append(k)
	    
	    # Build chr physical map by combining wild and evolved ones.
	    diploid_cells[n_d_c] = build_chromosome_physical_map(evolved_cells[sc],PM_hap,PM[sc])
	    
	    # calculate recombination rate	
	    diploid_cells[n_d_c] = calculate_recombination_rate(diploid_cells[n_d_c],cMorgan_values)
	    
	    # Calculating the Cell divison time and re-pack the cells for mitotic division.
	    # GD_hap_new.append([]);PM_hap_new.append([])
	    
	    genotypes,n_cell,PM_new,PM_hap_new,PM_strand = calculate_cell_division_time_diploid(diploid_cells[n_d_c],PM_fitness,genotypes,n_cell,PM_new,PM_hap,PM_strand)
	    
	    n_d_c += 1
	    
	    # Get the GD and PM fitness for the cells selected
	    #GD_fitness_new  = get_fitness_from_ID(GD_new,GD_fitness,GD_fitness_new,n_cell)
	    #PM_fitness_new  = get_fitness_from_ID(PM_new,PM_fitness,PM_fitness_new,n_cell)
	    k += 1
	    # n_cell += 1
    #for i1,i2 in enumerate(PM_new):
    #	print "PLACE II:",i1,"\t",i2
    
    return genotypes,n_cell,PM_new,PM_fitness,PM_hap,s_cells,PM_strand

def calculate_mating_params(evolved_cells,p_m):
    n_mate = (len(evolved_cells)/2) * p_m
    n_ec = (len(evolved_cells)/2) * (1-p_m)
    n_wt = (len(evolved_cells)/2) * (1-p_m)
    
    tot = n_mate + n_ec + n_wt
    n1 = n_mate/tot
    n2 = n_ec/tot
    n3 = n_wt/tot
    
    n2 = n3 + n2 
    
    return n1,n2,n3,tot
    
###########################################################################################################
# Just a program to executing the copy commands
###########################################################################################################
def copy_sub(n_cell,genotypes,evolved_cells,PM_new,PM):
    # Copying the properties of the whole cells
    genotypes[n_cell] = evolved_cells
    genotypes[n_cell][0] = genotypes[n_cell][1]

    #GD_new.append(GD)
    PM_new.append(PM)
    
    return genotypes,PM_new

##########################################################################################################
# To calculate the frequencies of genetic variations.
##########################################################################################################
def calculate_frequency_2(types,text,n_repeat,n_backcross,GDs,IDs,ID_fitness,param_text,new_ids):
    print GDs
    print IDs
    print ID_fitness
    sys.exit()
    # Combining new and old ids.    
    for i in xrange(len(IDs)):
	n = 0;tot = 0; txt = "B4BC"
	for j in xrange(len(GDs)):
	    tot += 1
	    for k in xrange(len(GDs[j])):
		if IDs[i] == int(GDs[j][k]):
		    n += 1
	
	# Finding the fitness type
	id = IDs[i]
	if ID_fitness[id] > 0:
	    fitness_type = "Positive"
	elif ID_fitness[id] < 0:
	    fitness_type = "Negative"
        else:
	    fitness_type = "Neutral"
	m = 0
	
	for l in xrange(len(new_ids)):
	    if IDs[i] == new_ids[l]:
		m = 1; txt = "DBC"
		
	if n == 0 or tot == 0:    
	    print text,"\t",n_repeat,"\t",n_backcross,'\t',types,"\t",IDs[i],"\t",fitness_type,"\t","0.0","\t",ID_fitness[id],"\t",param_text,"\t",txt 
	else:
	    print text,"\t",n_repeat,"\t",n_backcross,'\t',types,"\t",IDs[i],"\t",fitness_type,"\t",(n/tot),"\t",ID_fitness[id],"\t",param_text,"\t",txt

    return 0

##########################################################################################################
# To calculate the frequencies of genetic variations.
##########################################################################################################
def calculate_frequency(types,text,n_repeat,n_backcross,GDs,IDs,GD_fitness,param_text,s_size):
    """
    for i in IDs:
	print "IDs:",i
	
    for i1,i2 in GDs.iteritems():
	print "GDs:",i1,i2
	
    print "IDs:",len(IDs),IDs[0]
    print "GDs:",len(GDs),GDs[1]
    sys.exit()
    """
    # Combining new and old ids.
    for i in IDs:
	n = 0;tot = 0; txt = "B4BC"
	for j in xrange(s_size):
	    tot += 1
	    if i in GDs[j]:
		n += 1

	# Finding the fitness type
	if GD_fitness[i] > 0:
	    fitness_type = "Positive"
	elif GD_fitness[i] < 0:
	    fitness_type = "Negative"
        else:
	    fitness_type = "Neutral"
	
	if n == 0 or tot == 0:    
	    print text,"\t",n_repeat,"\t",n_backcross,'\t',types,"\t",i,"\t",fitness_type,"\t","0.0","\t",GD_fitness[i],"\t",param_text,"\t",txt 
	else:
	    print text,"\t",n_repeat,"\t",n_backcross,'\t',types,"\t",i,"\t",fitness_type,"\t",(n/tot),"\t",GD_fitness[i],"\t",param_text,"\t",txt

    return 0

##########################################################################################################
# To calculate the frequencies of genetic variations.
##########################################################################################################
#def calculate_frequency_1(types,text,n_repeat,n_backcross,GDs,IDs,ID_fitness,param_text,new_ids,s_cells)
def calculate_frequency_1(types,text,n_repeat,n_backcross,GDs,IDs,ID_fitness,param_text,new_ids,s_cells):
    
    # Adding only the selected cells
    GDs_new = []
    for i in s_cells:
	GDs_new.append(GDs[i])
	
    # Combining new and old ids.    
    IDs  = combine_markers(GDs_new,IDs)
        
    for i in IDs:
	n = 0;tot = 0; txt = "B4BC"
	for j in GDs_new:
	    tot = tot + 1
	    for k in j:
		if i == int(k):
		    n = n + 1
	
	# Finding the fitness type
	
	if ID_fitness[i] > 0:
	    fitness_type = "Positive"
	elif ID_fitness[i] < 0:
	    fitness_type = "Negative"
        else:
	    fitness_type = "Neutral"
	m = 0
	
	if i in new_ids:
	    m = 1; txt = "DBC"
	"""    
	for l in range(len(new_ids)):
	    if IDs[i] == new_ids[l]:
		m = 1; txt = "DBC"
	"""
	
	if n == 0 or tot == 0:    
	    print text,"\t",n_repeat,"\t",n_backcross,'\t',types,"\t",i,"\t",fitness_type,"\t","0.0","\t",ID_fitness[i],"\t",param_text,"\t",txt 
	else:
	    print text,"\t",n_repeat,"\t",n_backcross,'\t',types,"\t",i,"\t",fitness_type,"\t",(n/tot),"\t",ID_fitness[i],"\t",param_text,"\t",txt
    	    
    return 0

##########################################################################################################
# To Combine the new markers with OLD markers
##########################################################################################################
def combine_markers(new,old):
    new_temp = []
    # Combining the new ones
    for i in range(len(new)):
	for j in range(len(new[i])):   
	    new_temp.append(new[i][j])

    # Comibing the new ones with old ones
    for i in range(len(old)):
	new_temp.append(old[i])
    
    unique_ids = unique_grep(new_temp)
    
    return unique_ids

###############################################################################################
# To extract the unique mutations from the given data
###############################################################################################    
def unique_grep(data):
    set = {} 
    map(set.__setitem__, data, []) 
    
    return set.keys()
    
##########################################################################################################
# To catalogue the fitness type to track those alleles.
##########################################################################################################
def find_fitness_type(GD,GD_fitness):
    print GD
    print GD_fitness
    sys.exit()
    fitness_type = []
    for i in xrange(len(GD_fitness)):
	fitness_type.append([])
	for j in range(len(GD_fitness[i])):
	    if GD_fitness[i][j] == 0:
		type = "Neutral"
	    elif GD_fitness[i][j] < 0:
		type = "Negative"
	    else:
		type = "Positive"
		
	    fitness_type[i].append(type)
    
    return fitness_type
    
##########################################################################################################
# To get the fitness associated with genetic variations.
##########################################################################################################
def get_fitness_from_ID(GD_new,GD_fitness,GD_fitness_new,n_cell):
    GD_fitness_new.append([])
    
    for i in xrange(len(GD_new[n_cell-1])):
	id = int(GD_new[n_cell-1][i])
	if GD_fitness.has_key(id):
	    GD_fitness_new[n_cell-1].append(GD_fitness[id])

    return GD_fitness_new	    

##########################################################################################################
# Choose the viable cells after sporulation. And the best cell will be chosen in the same function
##########################################################################################################
def cell_viability(genotypes,n_cell,GD,GD_fitness,GD_hap,PM,PM_fitness,PM_hap):
    import random
    
    GD_new = [];GD_fitness_new = [];GD_hap_new = []
    PM_new = [];PM_fitness_new = [];PM_hap_new = []
    # Cells viability is 80-90% after sporulation
    n_viable_cell  = int(round((random.randrange(80,90)/100) * n_cell))
    n_delete = n_cell - n_viable_cell

    best_cell_id = 0;CDT = genotypes[0][1] 
    # Choosing the viable cells as well as the BEST CELL among the viable cells.
    for i in xrange(n_viable_cell):
	id = random.randrange(0,n_cell)
	if CDT > genotypes[id][1]:
	    best_cell_id = id
	    CDT = genotypes[id][1]

    # Make the best cell
    n_best_cell_copy = random.randrange(20,25)
    selected_cells = ones((n_best_cell_copy,5),dtype=int32)
         
    # Make the Best cell duplicate and make them ready for next round of mating
    for i in xrange(n_best_cell_copy):	
	selected_cells[i] = genotypes[best_cell_id]     
	# Copying the Gene duplication information
	GD_new.append(GD[best_cell_id])
	GD_fitness_new.append(GD_fitness[best_cell_id])
	GD_hap_new.append(GD_hap[best_cell_id])
	
	# Copying the point mutation information
	PM_new.append(PM[best_cell_id])
	PM_fitness_new.append(PM_fitness[best_cell_id])
	PM_hap_new.append(PM_hap[best_cell_id])    
    
    return selected_cells,GD_new,GD_fitness_new,GD_hap_new,PM_new,PM_fitness_new,PM_hap_new

##########################################################################################################
# Choose the viable cells after sporulation. And CHOOSE ABOUT 20-25 cells instead of ONLY ONE CELL.
##########################################################################################################
def cell_viability_1(genotypes,n_cell,GD,GD_fitness,GD_hap,PM,PM_fitness,PM_hap):
    import random
    
    GD_new = [];GD_fitness_new = [];GD_hap_new = []
    PM_new = [];PM_fitness_new = [];PM_hap_new = []
    # Cells viability is 80-90% after sporulation
    n_viable_cell  = int(round((random.randrange(80,90)/100) * n_cell))
    n_delete = n_cell - n_viable_cell

    best_cell_id = 0;CDT = genotypes[0][1]
    
    # Number cells to be chosen
    n_best_cell_copy = random.randrange(20,25)
    selected_cells = ones((n_best_cell_copy,5),dtype=int32)
    
    # Choosing the viable cells as well as the 20-25 cells for the next roud of mating among the viable cells.
    for i in xrange(n_best_cell_copy):
	id = random.randrange(0,n_cell)
	selected_cells[i] = genotypes[id]     
	# Copying the Gene duplication information
	GD_new.append(GD[id])
	GD_fitness_new.append(GD_fitness[id])
	GD_hap_new.append(GD_hap[id])
	
	# Copying the point mutation information
	PM_new.append(PM[id])
	PM_fitness_new.append(PM_fitness[id])
	PM_hap_new.append(PM_hap[id])
    
    return selected_cells,GD_new,GD_fitness_new,GD_hap_new,PM_new,PM_fitness_new,PM_hap_new

##########################################################################################################
# Choose the viable cells after sporulation. And It will return all the viable cells.
##########################################################################################################
def cell_viability_2(genotypes,n_cell,GD,GD_fitness,GD_hap,PM,PM_fitness,PM_hap,min_sample_size,max_sample_size):
    import random
    
    GD_new = [];GD_fitness_new = [];GD_hap_new = []
    PM_new = [];PM_fitness_new = [];PM_hap_new = []
    # Cells viability is 80-90% after sporulation
    n_viable_cell  = int(round((random.randrange(80,90)/100) * n_cell))

    best_cell_id = 0;CDT = genotypes[0][1]
    
    # Number cells to be chosen
    n_best_cell_copy = random.randrange(min_sample_size,max_sample_size)
    selected_cells = ones((n_viable_cell,5),dtype=int32)
    
    # Choosing the viable cells as well as the 20-25 cells for the next round of mating among the viable cells.
    for i in xrange(n_viable_cell):
	id = random.randrange(0,n_cell)
	selected_cells[i] = genotypes[id]     
	# Copying the Gene duplication information
	GD_new.append(GD[id])
	GD_fitness_new.append(GD_fitness[id])
	GD_hap_new.append(GD_hap[id])
	
	# Copying the point mutation information
	PM_new.append(PM[id])
	PM_fitness_new.append(PM_fitness[id])
	PM_hap_new.append(PM_hap[id])
    
    return selected_cells,GD_new,GD_fitness_new,GD_hap_new,PM_new,PM_fitness_new,PM_hap_new

##########################################################################################################
# Combine GD ids and their fitness effects
##########################################################################################################
def combine_fitness_2(GD,fitness):
    """
    for i in range(len(GD)):
	if len(GD[i]) != len(fitness[i]):
	    print "***",GD[i],fitness[i]
    """

    combined = {}
    for i in xrange(len(GD)):
	for j in xrange(len(GD[i])):
	    id = int(GD[i][j]);
	    value = fitness[i][j]
	    if combined.has_key(id):
	       continue
	    else:
	       combined[id] = value
    
    return combined

##########################################################################################################
# Combine GD ids and their fitness effects
##########################################################################################################
def combine_fitness(GD,fitness):    
    combined = {}
    k1 = 0; k2 = 0
    for i in GD:
	k1 += 1; k2 = 0	
	for j in i:
	    k2 += 1
	    id = int(j)
	    value = fitness[k1-1][k2-1]
	    if combined.has_key(id):
	       continue
	    else:
	       combined[id] = value
    
    return combined

##########################################################################################################
# Combine GD ids and their fitness effects with the OLD ones
##########################################################################################################
def combine_fitness_1(whole_fitness,GD,fitness):
    new_ids = []
    for i in GD.itervalues():
	for j in i:
	    ID = int(j)
	    if whole_fitness.has_key(ID):
	       pass
	    else:
	       whole_fitness[ID] = fitness[ID]
	       new_ids.append(ID)
	
    return whole_fitness, new_ids

##########################################################################################################
# this will build the physical map of chromosomes with mutations
##########################################################################################################
def build_chromosome_physical_map_1(evolved_cell,GD,PM,GD_id,PM_id):
    GD_hap = [];PM_hap = []
    GD_hap = GD; PM_hap = PM

    chr = array([0],dtype=float64)

    # NM Number of markers you want.
    NM = 10

    # Yeast Genome Phy Map
    chr_phy_map = array([resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
			 resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM)])
    

    # Map the GD Haplotype
    if len(GD_hap) == 0:
	le = 0
    elif len(GD_hap) < 2:
	# Get the Chr Number and Position
	chr_no,chr_position = split_hap(GD_hap[0])
	# 0.1 will indicate that its GD
	chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,GD_id[0],0.1)
    else:
	for i in range(len(GD_hap)):
	    ## Get the Chr Number and Position
	    chr_no,chr_position = split_hap(GD_hap[i])
	    # 0.1 will indicate that its GD
	    chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,GD_id[i],0.1)
    
    # Map the PM Haplotype
    if len(PM_hap) == 0:
	le = 0
    elif len(PM_hap) < 2:
	# Get the Chr Number and Position
	chr_no,chr_position = split_hap(PM_hap[0])
	# 0.2 will indicate thatt its PM
	chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,PM_id[0],0.2)
    else:
	for i in range(len(PM_hap)):
	    # Get the Chr Number and Position
	    chr_no,chr_position = split_hap(PM_hap[i])
	    # 0.2 will indicate thatt its PM
	    chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,PM_id[i],0.2)    

    return chr_phy_map

##########################################################################################################
# this will build the physical map of chromosomes with mutations
# Here the MAP Is STRING - This is the diff between the previous function which does the same.
##########################################################################################################
def build_chromosome_physical_map(evolved_cell,PM_hap,PM_id):
    
    
    #chr = array([0],dtype='|S16')
    chr = array([0,0,0,0,0,0,0,0,0,0],dtype='|S16')
    chr_phy_map = array([chr,chr,chr,chr,
                            chr,chr,chr,chr,
                            chr,chr,chr,chr,
                            chr,chr,chr,chr,
          		    chr,chr,chr,chr,
                            chr,chr,chr,chr,
                            chr,chr,chr,chr,
                            chr,chr,chr,chr])    
    
    
    # Number of markers you want.
    NM = 10
    """
    # Yeast Genome Phy Map
    chr_phy_map = array([resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
			 resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM),
                         resize(chr,NM),resize(chr,NM),resize(chr,NM),resize(chr,NM)])
    #print chr_phy_map
    #sys.exit()
    
    # Map the GD Haplotype
    if len(GD_id) == 0:
	le = 0
    elif len(GD_id) < 2:
	# Get the Chr Number and Position
	#id = GD_id[0]
	chr_no,chr_position = split_hap(GD_hap[GD_id[0]])
	# 0.1 will indicate that its GD
	chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,GD_id[0],0.1)
    else:
	i1 = 0
	for i in GD_id:
	    ## Get the Chr Number and Position
	    chr_no,chr_position = split_hap(GD_hap[i])
	    # 0.1 will indicate that its GD
	    chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,GD_id[i1],0.1)
	    i1 += 1
    """
    # Map the PM Haplotype
    if len(PM_id) == 0:
	le = 0
    elif len(PM_id) < 2:
	#id = PM_id[0]
	# Get the Chr Number and Position
	chr_no,chr_position = split_hap(PM_hap[PM_id[0]])
	# 0.2 will indicate thatt its PM
	chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,PM_id[0],0.2)
    else:
	i1 = 0 
	for i in PM_id:
	    # Get the Chr Number and Position
	    chr_no,chr_position = split_hap(PM_hap[i])
	    # 0.2 will indicate thatt its PM
	    chr_phy_map = map_hap(chr_phy_map,chr_no,chr_position,PM_id[i1],0.2)
	    i1 += 1

    return chr_phy_map


##########################################################################################################
# Mapping the Chr number and position on genome
##########################################################################################################
def map_hap(chr_phy_map,chr_no,chr_position,id,indicator):
    # General format of chr position
    # First digit after "."  - Type of marker (PM or GD) 
    # Remaining digit are the ID of the marker.
    
    # Preparing the Chr position with indicators
    t = str(indicator) + str(id)
    chr_position = float(chr_position) +  float(t)
    s = len(str(id)) + 1
    d = "%."+ str(s) + "f"
    chr_position = d % chr_position
    
    # Finding the number of markers already on the Chr so that we can place the next marker
    NM = number_of_markers(chr_phy_map[chr_no*2-1])
    
    # chr_phy_map[chr_no*2-1][NM] = float(chr_position)
    chr_phy_map[chr_no*2-1][NM] = str(chr_position)
    chr_phy_map[chr_no*2-1].sort()
    
    return chr_phy_map
    
##########################################################################################################
# Number of Markers on the Chromosome already
##########################################################################################################    
def number_of_markers(chrs):
    NM = 0
    for i in chrs:
    	if i != '0':
    	    NM += 1
    
    #c = Counter(chrs)
    
    #return 10 - c['0']
    return NM
	    
##########################################################################################################
# Split the Haplotypes
##########################################################################################################
def split_hap(hap):
    """
    NORMAL PYTHON
    temp_list = re.split(":",str(hap))
    temp_var = re.split("{",temp_list[0])
    chr_no = int(temp_var[1])
    chr_position = re.split("}",temp_list[1])    
    """

    """
    Advanced Python
    """
    pattern = "\{(\d+)\:\s(\d+)"

    #print hap
    b = re.search(pattern, str(hap))
    #print b
    #sys.exit()
    return int(b.group(1)),int(b.group(2))    

##########################################################################################################
# This will calculate the physical distance between markers in the chr
# and using that it will calculate recombination fractions and do 
# recombination between homologous chromosomes.
##########################################################################################################
def calculate_recombination_rate(diploid_cell,cMorgan_values):
    physical_distance = []

    # Using the chrmosomes from evolved cells not from the WT
    chr_no = 0 
    for i in xrange(1,32,2):
	t = 0
	NM = number_of_markers(diploid_cell[i])
	chr_no = chr_no + 1
	if NM > 1:
	    # Get the markers to get the distance between them
	    markers = bring_markers(diploid_cell[i])
	    
 	    #print "####################################MARKERS"
	    #print markers
	    #print "####################################MARKERS"
	    markers = sort_markers(markers)
 	    #print "####################################MARKERS"
	    #print markers
	    #print "####################################MARKERS"	    
	    #sys.exit()
            # Distance (physical) between the markers
	    physical_distance = find_distance(markers,chr_no,cMorgan_values)
	    
	    # Packing the Homolog Chromosomes into one pack to send it to recombine 
	    homolog_chr = chr_packing(diploid_cell[i-1],diploid_cell[i])
	    
            #a = array([0,0],dtype='|S16')
	    #temp_chr = array([resize(chr,(10,2))])
	    

	    # Calling the function to recombine
	    #print "BEFORE:",homolog_chr
	    #print physical_distance
	    temp_chr = meiosis.recombine(homolog_chr,physical_distance)
	    #print "AFTER:",temp_chr
	    #sys.exit()
	    
            # Re-Mapping the markers after recombination
	    diploid_cell[i-1],diploid_cell[i] = remap_markers(temp_chr,diploid_cell[i-1],diploid_cell[i])
    
    return diploid_cell

##########################################################################################################
# Pack the Homolog Chromosome into a desired format
##########################################################################################################
def chr_packing(chr_1,chr_2):

    chr = array([0,0],dtype='|S16')
    temp = array([resize(chr,(10,2))])
    # Funny thing to do. But it works in this way
    chr_map = temp[0]
    
    # Mapping the markers
    for i in xrange(len(chr_1)):
	chr_map[i][1]= chr_2[i]

    return chr_map

def sort_markers(chr_2):
    non_zeros = []
    for i in chr_2:
	if float(i) !=0:
	    non_zeros.append(float(i))
    
    non_zeros.sort()
    NZ = 0
    for i in xrange(len(chr_2)):
	if float(chr_2[i]) != 0:
	    chr_2[i] = str(non_zeros[NZ])
	    NZ += 1
    
    return chr_2
    
##########################################################################################################
# Pack the Homolog Chromosome into a desired format
##########################################################################################################
def remap_markers(temp_chr,chr_1,chr_2):

    for k in xrange(len(chr_1)):
	chr_1[k] = str(temp_chr[k][0])
	chr_2[k] = str(temp_chr[k][1])

    return chr_1,chr_2
    
##########################################################################################################
# Get the markers on the given Chromosome
##########################################################################################################
def bring_markers(chr):
    markers = [i for i in chr if i != 0]
    
    """
    markers = [] 
    for i in xrange(len(chr)):
	if chr[i] != 0:
	    markers.append(chr[i])
    """
    
    return markers 

##########################################################################################################
# Distance between the markers and convert it into cM/Kbps
##########################################################################################################
def find_distance(markers,chr_no,cMorgan_values):
    physical_distance = []

    for j in xrange(len(markers)-1):
	#print int(float(markers[j+1])), int(float(markers[j]))
	distance =  int(float(markers[j+1])) - int(float(markers[j]))
	
	# Finding the physical distance and dividing by 1000 to convert it to cM/Kbps
	# and multiplying by the cM value
	temp = (distance/1000) * cMorgan_values[chr_no-1]
	physical_distance.append(temp)
	
	# Calling the invHaldane function
	physical_distance[j] = meiosis.invHaldane(physical_distance[j])

    return physical_distance

##########################################################################################################
# Here the DNA replicates and the cells divide twice to make 4 haploids.
##########################################################################################################
def dna_replication_and_division_one(diploid_cell,PM_fitness,genotypes,n_cell,PM_new,PM_hap_new):

    chr = array([0,0],dtype='|S16')
    haploid_1 = resize(chr,(16,10))
    haploid_2 = resize(chr,(16,10))

    """
    cell division-I
    Each chromosome is replicated. So, for each chromosome
    We should choose the copy randomly. It will produce two cells
    Including the replication copy, there are four copies of each chr strand.
    We have to group them as two pairs randomly.
    1-wildtype original, 2-evolved original, 3-wiltype copy,4-evolved copy
    randomly choosing chromosome pairs
    """

    for i in xrange(0,2):
	# Deciding the Chr numbers for Haploids
	for j in range(16):
	    # Choose the strand
	    strand = random.randint(0,2)
	
	    # Decide the Chrosomes for each haploids
	    chr_no_1,chr_no_2 = decide_strands(j,strand)
	    
	    haploid_1[j] = diploid_cell[chr_no_1]	    
	    haploid_2[j] = diploid_cell[chr_no_2]
        
	# calculate the cell division time and pack the cell to divide again.

	# Appending new array to insert the Hap of new Haplotype.
	PM_hap_new.append([])
	genotypes,n_cell,PM_new,PM_hap_new = calculate_cell_division_time_haploid(haploid_1,PM_fitness,genotypes,n_cell,PM_new,PM_hap_new)

	# Appending new array to insert the Hap of new Haplotype.
	PM_hap_new.append([])
	genotypes,n_cell,PM_new,PM_hap_new = calculate_cell_division_time_haploid(haploid_2,PM_fitness,genotypes,n_cell,PM_new,PM_hap_new)	

	
    return genotypes,n_cell,GD_new,PM_new,GD_hap_new,PM_hap_new

##########################################################################################################
# Decide the Chromosome strands for haploids
##########################################################################################################
def decide_strands(chr,strand):
    chr = chr + 1

    chr_no_1 = chr*2-1-strand

    if strand == 0:
       strand = 1
    else:
        strand = 0

    chr_no_2 = chr*2-1-strand
         
    return chr_no_1,chr_no_2

##########################################################################################################
# Calculate cell division time and pack the cells for Mitotic dividing
##########################################################################################################
def calculate_cell_division_time_haploid(haploid,PM_fitness,genotypes,n_cell,PM_new,PM_hap_new):

    max_CDT = 180
    n_chr = 16
    # Get the list of marker ids in the Haploid cell
    PM,PM_hap_new = split_markers_id(haploid,PM_hap_new,n_cell,n_chr)

    # Adding them to new storage place
    #GD_new.append(GD)
    PM_new.append(PM)
        
    # Calculate total fitness effects due to markers
    fitness = total_fitness(PM,PM_fitness)

    # Cell division time
    CDT = cell_division_time(max_CDT,fitness)

    # Build the new cell to divide mitotically
    genotypes,n_cell = cell_building(genotypes,CDT,PM,n_cell)
    
    return genotypes,n_cell,PM_new,PM_hap_new


##########################################################################################################
# Calculate cell division time and pack the cells for Mitotic dividing
##########################################################################################################
def calculate_cell_division_time_diploid(diploid,PM_fitness,genotypes,n_cell,PM_new,PM_hap_new,PM_strand):
    
    max_CDT = 180
    n_chr = 32
    
    # Get the list of marker ids in the Haploid cell
    PM,PM_hap_new,PM_strand_new = split_markers_id(diploid,PM_hap_new,n_cell,n_chr)

    # Adding them to new storage place
    #GD_new.append(GD)
    PM_new.append(PM)
    
    # Adding the strand information
    # GD_strand.append(GD_strand_new)
    # PM_strand.append(PM_strand_new)
    
    # Calculate total fitness effects due to markers
    fitness = total_fitness(PM,PM_fitness)
    
    # Cell division time
    CDT = cell_division_time(max_CDT,fitness)
    
    # Build the new cell to divide mitotically
    genotypes,n_cell = cell_building(genotypes,CDT,PM,n_cell)
    
    return genotypes,n_cell,PM_new,PM_hap_new,PM_strand

##########################################################################################################
# To convert the matrix of floats to strings to avoid precision changing problem.
# Background: when you send a float to another funtion it automatically converted to some other precision
# To avoid that I have to convert it to string before I send.
########################################################################################################## 
def convert_to_matrix_string(haploid):
    chr = array([0,0],dtype='|S16')
    haploid_1 = resize(chr,(16,10))
    for i in xrange(16):
	for j in xrange(len(haploid[i])):
	    haploid_1[i][j] = str(haploid[i][j])
    
    return haploid_1

##########################################################################################################
# To calculate the Cell division time based on fitness effect
##########################################################################################################     
def cell_building(genotypes,CDT,PM,n_cell):
    # Assigning next cell division time
    genotypes[n_cell][0] = CDT
    # Assigning current cell division time of th cells.
    genotypes[n_cell][1] = CDT
    # Number of possible cell divisions
    genotypes[n_cell][2] = 16
    # Number of point mutations
    genotypes[n_cell][3] = len(PM)
    #Number of Gene duplications
    genotypes[n_cell][4] = 0 # len(GD)
    # Adding the type of the cells(1 - haploid, 2 - diploid)
    genotypes[n_cell][5] = 2
    
    # Number of cells
    n_cell += 1
    
    
    return genotypes,n_cell
    
##########################################################################################################
# To calculate the Cell division time based on fitness effect
########################################################################################################## 
def cell_division_time(cell_division_time,alpha):
    cell_division_max = 180
    cell_division_min = 90
    
    # Cell Division time to "Alpha"-refer the coding documents.
    beta = cell_division_max - cell_division_min
    natural_log = -(math.log(2) * (alpha))
    CDT = math.ceil(beta * math.exp(natural_log) + cell_division_min)
    
    return CDT    

##########################################################################################################
# To calculate the total fitness from the markers
########################################################################################################## 
def total_fitness(PM,PM_fitness):
    fitness = 0
    
    # Total Fitness from GD
    
    #for i in GD:
    #	fitness += GD_fitness[int(i)]
    
    for i in PM:
	fitness += PM_fitness[int(i)]
    
    return fitness
##########################################################################################################
# Trick ? or Hack ? have to find a solution for this.
##########################################################################################################
def find_correct_id(id,fitness):
    n = 0 
    while n == 0:
      id = id * 10
      if fitness.has_key(id):
	n = 1

    return id   

##########################################################################################################
# Extract the GD, PMs ids from Chr map
##########################################################################################################    
def split_markers_id(haploid_1,PM_hap_new,n_cell,n_chr):

    PM = []; PM_strand = {}
    
    for i in xrange(n_chr):
	for j in haploid_1[i]:
	    if float(j) != 0:
	       # Get Type and ID
	       geno_type,geno_id = split_id_type(j)
	       
	       # Store the chr type for PM and GD
	       PM_strand= store_chr_strand(geno_type,geno_id,PM_strand,i)
	       
	       # Build GD_hap and PM_hap
	       # GD_hap_new, PM_hap_new = build_hap(haploid_1[i][j],i,GD_hap_new,PM_hap_new,n_cell)
	       
	       # Store the markers in right category
	       PM = store_markers(PM,geno_type,geno_id)

    return PM,PM_hap_new,PM_strand


##########################################################################################################
# Extract the GD, PMs ids from Chr map
########################################################################################################## 
def store_chr_strand(geno_type,geno_id,PM_strand,i):
    # Unary Operator
    #strand = 2 if i % 2 == 0 else 1
    strand = 2
    
    if PM_strand.has_key(geno_id):
	pass
    else:	    
	PM_strand[geno_id] = strand
	
    """
    if geno_type == 1:
	if GD_strand.has_key(geno_id):
	    pass
	else:	    
	   GD_strand[geno_id] = strand
    else:
	if PM_strand.has_key(geno_id):
	    pass
	else:	    
	   PM_strand[geno_id] = strand
    """
    
    return PM_strand

##########################################################################################################
# Extract the GD, PMs ids from Chr map
##########################################################################################################    
def split_id_type(marker):
    pattern = "(\d+)\.(\d)(\d+)"
    b = re.search(pattern, str(marker))

    return int(b.group(2)),int(b.group(3))

##########################################################################################################
# Extract the GD, PMs ids from Chr map and BUILD GD_HAP and PM_HAP
##########################################################################################################     
def build_hap(marker,chr_n,PM_hap_new,n_cell):
    
    pattern = "(\d+)\.(\d)(\d+)"
    b = re.search(pattern, str(marker))
    
    hap = "{" + str(chr_n+1)+ ": " +str(b.group(3)) +"}"
    
    # Check whether its GD and PM and append    
    if int(b.group(2)) == 1:
	GD_hap_new[n_cell].append(str(hap))
    else:
	PM_hap_new[n_cell].append(str(hap))

    print GD_hap_new, PM_hap_new
    sys.exit()
    return GD_hap_new, PM_hap_new
    
##########################################################################################################
# Store them in to the right category
##########################################################################################################     
def store_markers(PM,type,id):
    PM.append(id)   
	
    return PM

##########################################################################################################
# Convert the data form from List of lists to Dict of lists
##########################################################################################################
def convert_format(data):
    GD = {}
    k = 0 
    for i in data:
	GD[k] = [int(j) for j in i]
	k += 1
	
    return GD

##########################################################################################################
# Convert the data form from Dict of lists to List of lists 
##########################################################################################################
def convert_format_1(data):
    
    GD = [ data[i] for i in data]
    """  
    for i in xrange(len(GD)):
	t = [ fitness[j] for j in GD[i]]
	fitness_temp.append(t)
    """
    return GD

##########################################################################################################
# 
##########################################################################################################
def convert_format_2(PM,PM_hap):
    hap = {}
    for i in PM:
	#temp = [PM_hap[j].strip() for j in i]
	for j in i:
	    hap[j] = PM_hap[j-1].strip()
    
    return hap
	
##########################################################################################################
# Expand the data
##########################################################################################################
def expand_format(evolved_cells,n_generations,N_BN,n_cell):
    size =  n_cell * 2**n_generations 
    genotypes = ones((size,6),dtype=int32)
    for i in xrange(int(n_cell)):
	genotypes[i] = evolved_cells[i]
	
    return genotypes

##########################################################################################################
# Check for the types of Epistasis
##########################################################################################################
def check_epistasis_type(CDT,genotypes,n_cell,epistasis,j,BC):
    start = n_cell-3
    
    for i in xrange(start-1,n_cell):
	if genotypes[i][1] < CDT:
	    epistasis.write('%d\t%d\t%d\t%s\n' % (j, BC, i, "Positive Epistasis"))
	else:
	    epistasis.write('%d\t%d\t%d\t%s\n' % (j, BC, i, "Negative Epistasis"))
	    
    return 0

##########################################################################################################
# Checking the PM fitness and IDs to fix a bug 
##########################################################################################################	
def check_function(PM_2_BC,PM_fitness_2_BC):
    
    #print PM_2_BC
    #print PM_fitness_2_BC
    
    #print "INSIDE CHECK_FUNCTION"
    for i1,i2 in PM_2_BC.iteritems():
	for j in i2:
	    if PM_fitness_2_BC.has_key(j):
		pass 
	    else:
		print "####:",j
		print "#################################"
		print PM_2_BC
		print "#################################"
		print PM_fitness_2_BC
		print "#################################"
		sys.exit()
	    
    """
    print "###########################################"
    print PM_2_BC
    print "###########################################"
    print PM_fitness_2_BC
    print "###########################################"
    print PM_hap_2_BC
    print "###########################################"
    """

def check_function_1(PM_2_BC,PM_fitness_2_BC):    

    #print "INSIDE CHECK_FUNCTION"
    for i in PM_2_BC:
	for j in i:
	    if PM_fitness_2_BC.has_key(j):
		pass 
	    else:
		print "Missing Key:",j
		
    return 0

##########################################################################################################
# Converts diploid cells to Haploid cells 
##########################################################################################################
def diploids_2_haploids(evolved_cells,PM,PM_fitness,PM_strand,s_size):
    
    # Getting the index of DIPLOID cells
    dip_index = get_diploid_cells_index(evolved_cells,s_size)

    # 99% of Diploids are forced into Haploids
    n_cells = int(round(len(dip_index) * 0.99))

    # n =  99% of no. of diploids * 2 + no. of haploids + 1% of no. of diploids
    n = (n_cells * 2) + (len(evolved_cells) - len(dip_index)) + (len(dip_index) - n_cells)
    
    PM_temp = [];
    genotypes = ones((n,6),dtype=int32)

    choosen_cells = random.sample(dip_index,n_cells)
    k = -1
    k1 = 0 
    for i in xrange(int(s_size)): #evolved_cells:
	# Force the selected Diploids to Haploids
	if k1 in choosen_cells: # and evolved_cells[i][5] == 2:
	    k += 2
	    (genotypes[k-1],genotypes[k],PM1_temp,PM2_temp) = build_haploid_cells(evolved_cells[i],PM[k1],PM_fitness,PM_strand)
	    
	    # Copying the GDs and PMs
	    #GD_temp.append(GD1_temp); 
	    PM_temp.append(PM1_temp); 
	    #GD_temp.append(GD2_temp); 
	    PM_temp.append(PM2_temp); 
	    
	# Remaining cells will remain as Diploids
	else:
	    k += 1
	    genotypes[k] = evolved_cells[k1]
	    genotypes[k][0] = genotypes[k][1]
	    #GD_temp.append(GD[k1]) 
	    PM_temp.append(PM[k1]) 
	    
	k1 += 1

    return genotypes,PM_temp,PM_fitness

def get_diploid_cells_index(cells,s_size):
    
    return [i for i in xrange(int(s_size)) if cells[i][5] == 2]

##########################################################################################################
# Converts diploid cells to Haploid cells 
##########################################################################################################
def diploids_2_haploids_1(evolved_cells,GD,GD_fitness,GD_strand,PM,PM_fitness,PM_strand,n_evolved_cells):

    # Getting the index of DIPLOID cells
    dip_index = [i for i in range(len(evolved_cells)) if evolved_cells[i][5] == 2]
    
    # 99% of Diploids are forced into Haploids
    n_cells = int(round(len(dip_index) * 0.99))
    
    # n =  99% of no. of diploids * 2 + no. of haploids + 1% of no. of diploids
    n = (n_cells * 2) + (len(evolved_cells) - len(dip_index)) + (len(dip_index) - n_cells)
    
    GD_temp = []; GD_fitness_temp = [];GD_fitness_t = {}
    PM_temp = []; PM_fitness_temp = [];PM_fitness_t = {}
    GD_temp_t = {}; PM_temp_t = {}
    
    genotypes = ones((n,6),dtype=int32)
    
    choosen_cells = random.sample(dip_index,n_cells)
    k = -1
    for i in range(len(evolved_cells)):
	# Force the selected Diploids to Haploids
	if i in choosen_cells and evolved_cells[i][5] == 2:
	    k += 2
	    GD1_temp = []; GD1_fitness_temp = []; GD2_temp = []; GD2_fitness_temp = []
	    PM1_temp = []; PM1_fitness_temp = []; PM2_temp = []; PM2_fitness_temp = []
	    
	    (genotypes[k-1],genotypes[k],GD1_temp,GD1_fitness_temp,GD2_temp,GD2_fitness_temp,PM1_temp,PM1_fitness_temp,PM2_temp,PM2_fitness_temp) = build_haploid_cells(evolved_cells[i],GD[i],GD_fitness[i],GD_strand,PM[i],PM_fitness[i],PM_strand)
    
	    # Copying the GDs and PMs
	    GD_temp_t[k-1] = GD1_temp
	    GD_temp_t[k] = GD2_temp
	    PM_temp_t[k-1] = PM1_temp
	    PM_temp_t[k] = PM2_temp		    

	    #GD_temp.append(GD1_temp); #GD_fitness_temp.append(GD1_fitness_temp)
	    #PM_temp.append(PM1_temp); #PM_fitness_temp.append(PM1_fitness_temp)
	    #GD_temp.append(GD2_temp); #GD_fitness_temp.append(GD2_fitness_temp)
	    #PM_temp.append(PM2_temp); #PM_fitness_temp.append(PM2_fitness_temp)
	    
	    GD_fitness_t = copy_fitness(GD_fitness_t,GD1_temp,GD1_fitness_temp)
	    GD_fitness_t = copy_fitness(GD_fitness_t,GD2_temp,GD2_fitness_temp)
	    PM_fitness_t = copy_fitness(PM_fitness_t,PM1_temp,PM1_fitness_temp)
	    PM_fitness_t = copy_fitness(PM_fitness_t,PM2_temp,PM2_fitness_temp)	    

	# Remaining cells will remian as Diploids
	else:
	    k += 1
	    genotypes[k] = evolved_cells[i]
	    genotypes[k][0] = genotypes[k][1]
	    
	    GD_temp_t[k] = GD[i]
	    PM_temp_t[k] = PM[i]	    
	    #GD_temp.append(GD[i]); #GD_fitness_temp.append(GD_fitness[i])
	    #PM_temp.append(PM[i]); #PM_fitness_temp.append(PM_fitness[i])
	    
	    GD_fitness_t = copy_fitness(GD_fitness_t,GD[i],GD_fitness[i])
	    PM_fitness_t = copy_fitness(PM_fitness_t,PM[i],PM_fitness[i])
	    
    # STEP NEXT to Choose only 500 Cells
    tot_cells = n_cells * 2 + len(dip_index) - n_cells
    mm = random.sample(range(0,tot_cells),n_evolved_cells)
    n = n_evolved_cells * 2**5
    genotypes_s = ones((n,6),dtype=int32)
    GD_s = {}; PM_s = {}    
    for i in xrange(len(mm)):
	genotypes_s[i] = genotypes[mm[i]]
	GD_s[i] = GD_temp_t[mm[i]]
	PM_s[i] = PM_temp_t[mm[i]] 

    return genotypes_s,GD_s,GD_fitness_t,PM_s,PM_fitness_t


def copy_fitness(GD_fitness,GD,GD1_fitness):
    for i in xrange(len(GD)):
	if GD_fitness.has_key(GD[i]):
	    pass
	else:
	    GD_fitness[GD[i]] = GD1_fitness[i] 
    
    return GD_fitness

##########################################################################################################
# 
##########################################################################################################
def build_haploid_cells(evolved_cell,PM,PM_fitness,PM_strand):
    """
    k = 0 
    for i in PM_strand:
	print PM_strand[i]
	if PM_strand[i] == 2:
	    k += 1
    print "***",k
    sys.exit()
    """
    hap1 = [180,180,16,0,0,1]
    hap2 = [180,180,16,0,0,1]
    
    # Fitness for Hap2 and Calculate NEW CDT from fitness
    hap1,PM1_temp = total_fitness_haploid(PM,PM_fitness,PM_strand,1,hap1)

    # Fitness for Hap2 and Calculate NEW CDT from fitness
    hap2,PM2_temp = total_fitness_haploid(PM,PM_fitness,PM_strand,2,hap2)
    
    return hap1,hap2,PM1_temp,PM2_temp

##########################################################################################################
# Fitness calculation during DIPLOID to HAPLOID
##########################################################################################################
def total_fitness_haploid(GD,GD_fitness,GD_strand,ploidy,hap_cell):
    fitness = 0;n_gd = 0
    GD_temp = [] 
    for i in GD:
	if GD_strand[i] == ploidy:
	    # hap_cell[0] = calculate_cell_division_time(hap_cell[0],GD_fitness[i])
	    GD_temp.append(i)
	    n_gd += 1
	    
    # Update the next cell division time
    # hap_cell[1] = hap_cell[0]
    # Update the number of Mutations
    hap_cell[3] = n_gd
    
    return hap_cell,GD_temp

##########################################################################################################
# Special sampling Function for Haploid and Diploids mixed colonies
##########################################################################################################
def sampling_cells(genotypes,PM,PM_fitness,n_cell):
    n = n_cell
    #GD1 = {}; #GD_fitness1 = {}
    PM1 = {}; #PM_fitness1 = {}
    choosen_cells = random.sample(xrange(int(n_cell)),int(n))
    n = n * 2**5
    genotypes_1 = ones((n,6),dtype=int32)
    j = 0 
    for i in choosen_cells:
	# Choose only the Haploids
	genotypes_1[j] = deepcopy(genotypes[i])
	PM1[j] = PM[i]
	"""
	k1 = 0
	for k in GD1[j]:
	    GD_fitness1[k] = GD_fitness[i][k1]
	    k1 += 1
	    
	k1 = 0  
	for k in PM1[j]:
	    PM_fitness1[k] = PM_fitness[i][k1]
	    k1 += 1
	"""
	j += 1

    return genotypes_1,PM1,j
    
##########################################################################################################
# 
##########################################################################################################
def testing_baba(GD,fitness,k):
    for ii in range(len(GD)):
        for jj in range(len(GD[ii])):
            ID = GD[ii][jj]
            if fitness.has_key(ID):
                pass
            else:
                print "Place:",k,"MISSING KEY:",ID
    return 0

##########################################################################################################
# Make the cell division time as same for all the cells.
# Assumption: Cells divide more or less in same speed under normal conditions
##########################################################################################################
def prepare_for_non_s_exp(genotypes,n):
    for i in range(n):
	genotypes[i][0] = 180
	genotypes[i][1] = 180
    
    return genotypes

##########################################################################################################
# New Cell Division time is calculated based on the Fitness effects carried by PMs and GDs
##########################################################################################################
def prepare_for_s_exp(genotypes,PM,PM_fitness,PM_strand,n):
    #for j in range(n):
    #print genotypes[j]
    
    for i in xrange(n):
	fitness = 0
	
	# To add all the fitness effects
	# fitness += total_fitness_1(GD[i],GD_fitness)
	genotypes[i] = total_fitness_1(PM[i],PM_fitness,genotypes[i])
	genotypes[i][0] = genotypes[i][1]
	
	# Calculate the cell division time.
	# print i,"\t",genotypes[i],"\t",fitness
	# genotypes[i][1] = calculate_cell_division_time(genotypes[i][1],fitness)
	
    
    return genotypes

##########################################################################################################
# Summing up the fitness effects
##########################################################################################################
def total_fitness_1(IDs,fitness,genotypes):
    fitness_tot = 0
    for i in IDs:
	genotypes[1] = calculate_cell_division_time(genotypes[1],fitness[i])
	    
    return genotypes
    
    
def calculate_cell_division_time(current_CDT,alpha):
    # Cell Division time to "Alpha"-refer the coding documents.
    cell_division_max = 180
    cell_division_min = 90
    CDT_change = 0
    if alpha > 0:
        CDT_change = (current_CDT - cell_division_min)/(cell_division_max - cell_division_min)
        CDT_change *= alpha
        new_CDT = current_CDT - CDT_change
    else:
        new_CDT = current_CDT + abs(alpha)

    return  new_CDT

def freq_test(data,l):
    merged_data = [item for i in data for item in i]
    uniq_data = unique_grep(merged_data)
    freq = {}
    for i in uniq_data:
	k = 0
	for j in merged_data:
	    if i == j:
		k += 1
	freq[i] = k/l
	
    return freq

def freq_test1(data,l):
    merged_data = [item for i in data for item in data[i]]
    uniq_data = unique_grep(merged_data)
    freq = {}
    for i in uniq_data:
	k = 0
	for j in merged_data:
	    if i == j:
		k += 1
	freq[i] = k/l
	
    return freq

def freq_test2(data,uniq_data,l,fitness):
    merged_data = [item for i in data for item in i]
    #uniq_data = unique_grep(merged_data)
    freq = {}
    for i in uniq_data:
	k = 0
	for j in merged_data:
	    if i == j:
		k += 1
	freq[i] = k/l

    lost = 0;n_bene = 0
    for i in freq:
	if fitness[i] >0 and freq[i] == 0:
	    lost += 1
	if fitness[i] >0:
	    n_bene += 1
    
    lost_per = lost/n_bene
    
    return lost_per


def freq_test3(data,uniq_data,l,fitness):
    merged_data = [item for i in data for item in data[i]]
    #uniq_data = unique_grep(merged_data)
    freq = {}
    for i in uniq_data:
	k = 0
	for j in merged_data:
	    if i == j:
		k += 1
	freq[i] = k/l

    lost = 0;n_bene = 0
    for i in freq:
	if fitness[i] >0 and freq[i] == 0:
	    lost += 1
	if fitness[i] >0:
	    n_bene += 1
    
    lost_per = lost/n_bene
    
    return lost_per