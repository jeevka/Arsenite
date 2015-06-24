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
#import Artificial_Cells
###################################################################################################
# Global variables - Importing from main Modules.
###################################################################################################
# Gamma distribution values for Mutations. Refer the coding document.Or see the refernce Below.
# scale_alpha = 33
# shape_beta = 2
fitness_affecting_mutation_rate = Yeast_Simulator.fitness_affecting_mutation_rate
#shape = Yeast_Simulator.mutation_shape
#scale = Yeast_Simulator.mutation_scale
#mutation_fitness = Yeast_Simulator.mutation_fitness
#PM_haplotype = Yeast_Simulator.mutation_haplotype

###################################################################################################
# Introducing Mutations Based on Yeast Mutation rate.
###################################################################################################
def introduce_mutation(mother_cell,daughter_cell,beneficial_mutation_rate,k,tp,point_mutation,n_mut,PM,PM_fitness_Rate,PM_fitness_Lag,ploidy,fix_cal,CRT,CLT,N_HAP1,N_HAP2,N_HAP3,CH,N_CH,MID_Hap):
    
    # Decide whether its fitness altering mutation or neutral.
    # Beneficial Mutation rate is 575 out of 10000 fitness altering mutations.
    # Its 0.088 or 88/1000. Plz refer the coding value documents.
    # From Coding documents: 1.53% of all mutations are Fitness altering mutations.
    
    fitness_Rate = 0
    fitness_Lag = 0

    fitness_Rate = assign_haplotype()
        
    # Choose the cell:Parent or child to have mutation.
    ran_n = c_libc_random() % 0.99 + 0
    if ran_n <= 0.5:
        # Assigning the Chromosome and the position.
    
        # FPS1 LOF
        if fitness_Rate == 1:
            n_mut += 1
            N_HAP1 += 1
            # Calculate the RATE CDT
            if mother_cell[1] == 162:
                mother_cell[1] = 130 #Yeast_Simulator.calculate_cell_division_time(mother_cell[1],fitness_Rate)
                CRT[k] = mother_cell[1]
            
                # Calculate the Lag Time  
                CLT[k] = 276 #Yeast_Simulator.calculate_Lag_time(CLT[k],fitness_Lag)
                #CLT[k] = mother_cell[1]

            elif mother_cell[1] == 134:
                mother_cell[1] = 130 #Yeast_Simulator.calculate_cell_division_time(mother_cell[1],fitness_Rate)
                CRT[k] = mother_cell[1]
        
                # Calculate the Lag Time  
                CLT[k] = 276 #Yeast_Simulator.calculate_Lag_time(CLT[k],fitness_Lag)
                
            else:
                pass
            
            if CH.has_key(k):
                if CH[k] == "W":
                    CH[k] = "F"
                    N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
                else:
                    CH[k] = CH[k] + "F"
                    N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
            else:
                CH[k] = "F"
                N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
                
            mother_cell[5] = N_HAP1
            PM[k].append(n_mut)
            MID_Hap[n_mut] = "FPS1"        
        
        # ASK10 LOF 
        if fitness_Rate == 2:
            N_HAP2 += 1
            n_mut += 1
            if mother_cell[1] == 162:
                # Calculate the RATE CDT
                mother_cell[1] = 134 #Yeast_Simulator.calculate_cell_division_time(mother_cell[1],fitness_Rate)
                CRT[k] = mother_cell[1]
                
                # Calculate the Lag Time  
                CLT[k] = 503 #Yeast_Simulator.calculate_Lag_time(CLT[k],fitness_Lag)
                #CLT[k] = mother_cell[1]
            else:
                pass
                
            if CH.has_key(k):
                if CH[k] == "W":
                    CH[k] = "A"
                    N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
                else:
                    CH[k] = CH[k] + "A"
                    N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
            else:
                CH[k] = "A"
                N_CH[k] = N_HAP1 + N_HAP2
            
            mother_cell[4] = N_HAP2
            MID_Hap[n_mut] = "ASK10"
            PM[k].append(n_mut)


        # ACR3 Duplication 
        if fitness_Rate == 3:
            N_HAP3 += 1
            n_mut += 1
            if mother_cell[1] == 162 or mother_cell[1] == 134 or mother_cell[1] == 130:
                # Calculate the RATE CDT
                mother_cell[1] = 123 #Yeast_Simulator.calculate_cell_division_time(mother_cell[1],fitness_Rate)
                CRT[k] = mother_cell[1]
                
                # Calculate the Lag Time  
                CLT[k] = 630 #Yeast_Simulator.calculate_Lag_time(CLT[k],fitness_Lag)
                #CLT[k] = mother_cell[1]
            else:
                pass
                
            if CH.has_key(k):
                if CH[k] == "W":
                    CH[k] = "AC"
                    N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
                else:
                    CH[k] = CH[k] + "AC"
                    N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
            else:
                CH[k] = "AC"
                N_CH[k] = N_HAP1 + N_HAP2 + N_HAP3
            
            mother_cell[4] = N_HAP3
            MID_Hap[n_mut] = "ACR3"
            PM[k].append(n_mut)
            
        point_mutation.write("\n")
        mother_cell[3] += 1
        
        # Appending the mutation number
        # PM[k].append(n_mut)
        
    else:
        # FPS1 LOF
        if fitness_Rate == 1:
            n_mut += 1
            N_HAP1 += 1
            # Calculate the RATE CDT
            if daughter_cell[1] == 162:
                daughter_cell[1] = 130 #Yeast_Simulator.calculate_cell_division_time(daughter_cell[1],fitness_Rate)
                CRT[tp] = daughter_cell[1]
            
                # Calculate the Lag Time  
                CLT[tp] = 276 #Yeast_Simulator.calculate_Lag_time(CLT[tp],fitness_Lag)
                #CLT[tp] = daughter_cell[1]
                
            elif daughter_cell[1] == 134:
                daughter_cell[1] = 130 #Yeast_Simulator.calculate_cell_division_time(mother_cell[1],fitness_Rate)
                CRT[tp] = daughter_cell[1]
        
                # Calculate the Lag Time  
                CLT[tp] = 276 #Yeast_Simulator.calculate_Lag_time(CLT[k],fitness_Lag)                
            else:
                pass
                    
            if CH.has_key(tp):
                if CH[tp] == "W":
                    CH[tp] = "F"
                    N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
                else:
                    CH[tp] = CH[tp] + "F"
                    N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
            else:
                CH[tp] = "F"
                N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
            
            daughter_cell[5] = N_HAP1
            MID_Hap[n_mut] = "FPS1"
            PM[tp].append(n_mut)
        
        # ASK10 LOF
        if fitness_Rate == 2:
            N_HAP2 += 1
            n_mut += 1
            # Calculate the RATE CDT
            if daughter_cell[1] == 162:
                daughter_cell[1] = 134 #Yeast_Simulator.calculate_cell_division_time(daughter_cell[1],fitness_Rate)
                CRT[tp] = daughter_cell[1]
            
                # Calculate the Lag Time  
                CLT[tp] = 503 #Yeast_Simulator.calculate_Lag_time(CLT[tp],fitness_Lag)
                #CLT[tp] = daughter_cell[1]        
            else:
                pass

            if CH.has_key(tp):
                if CH[tp] == "W":
                    CH[tp] = "A"
                    N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
                else:
                    CH[tp] = CH[tp] + "A"
                    N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
            else:
                CH[tp] = "A"
                N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
                
            daughter_cell[4] = N_HAP2
            MID_Hap[n_mut] = "ASK10"
            PM[tp].append(n_mut)
        
        # ACR3 Duplication
        if fitness_Rate == 3:
            N_HAP3 += 1
            n_mut += 1
            # Calculate the RATE CDT
            if daughter_cell[1] == 162 or mother_cell[1] == 134 or mother_cell[1] == 130:
                daughter_cell[1] = 123 #Yeast_Simulator.calculate_cell_division_time(daughter_cell[1],fitness_Rate)
                CRT[tp] = daughter_cell[1]
            
                # Calculate the Lag Time  
                CLT[tp] = 630 #Yeast_Simulator.calculate_Lag_time(CLT[tp],fitness_Lag)
                #CLT[tp] = daughter_cell[1]        
            else:
                pass

            if CH.has_key(tp):
                if CH[tp] == "W":
                    CH[tp] = "AC"
                    N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
                else:
                    CH[tp] = CH[tp] + "AC"
                    N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
            else:
                CH[tp] = "AC"
                N_CH[tp] = N_HAP1 + N_HAP2 + N_HAP3
                
            daughter_cell[4] = N_HAP3
            MID_Hap[n_mut] = "ACR3"
            PM[tp].append(n_mut)

            
        point_mutation.write("\n")
        daughter_cell[3] +=  1
       
        # Appending the mutation number
        # PM[tp].append(n_mut)
    

    return mother_cell,daughter_cell,PM, CRT,CLT,N_HAP1,N_HAP2,N_HAP3,CH,N_CH,n_mut,MID_Hap
    
###################################################################################################
# Assigning the Haplotype structure.
###################################################################################################
def assign_haplotype():
        
    # Chromosome numbers and the length.
    #chrom = {
    #         1:230208,2:813178,3:316617,4:1531918,
    #         5:576869,6:270148,7:1090946,8:562643,
    #         9:439885,10:745745,11:666454,12:1078175,
    #         13:924429,14:784333,15:1091289,16:948062
    #         }
    #length = 0
    #for i in chrom:
    #    length += chrom[i]
    #print length
    #sys.exit()
    #haplotype_structure = {}
    # Choose a chromosome randomly.
    #chr_no = random.randint(1,16)
    # Choose a specific position in chromosome randomly.
    #bp_position = random.randint(1,chrom[chr_no])
    #haplotype_structure[chr_no] = bp_position
    
    hap = 0
    bps = random.randint(1,12070899)
    #bps = random.randint(1,120700)
    
    # Values for FPS1
    if bps <= 643:
        #print "FPS1"
        hap = 1

    # Values for ASK10
    elif bps > 643 and bps <= 1586:
        #print "ASK10"
        hap = 2
        
    # Values for ACR3
    elif bps > 1586 and bps <= 2529:
        #print "ACR10"
        hap = 3
    
        
    return hap
    
####################################################################################################    

"""
Reference:
Sarah B et al.  Spontaneous Mutations in Diploid Saccharomyces Cerevisiae: More Beneficial Than Expected. 2004 Genetics.
"""
####################################################################################################    
