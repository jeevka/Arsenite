from __future__ import division
import sys
import random
from copy import deepcopy
from heapq import *
#from heapq import heappush, heappop
from stdlib cimport *
from numpy import *
import numpy as np
cimport numpy as np
cimport cython
import gc
import re  
#from guppy import hpy

#dtype = np.int
#ctypedef np.int64_t dtype_t

import Yeast_Simulator
import Yeast_Simulator_Subprograms
import Yeast_Simulator_Mutations
import BackCrossing


#cdef int number_of_bottlenecks = Yeast_Simulator.number_of_bottlenecks
# cdef int n_individuals = Yeast_Simulator.n_individuals
cdef int cell_max_age = Yeast_Simulator.cell_max_age
cdef double mutation_rate = Yeast_Simulator.mutation_rate
cdef double gene_duplication_rate = Yeast_Simulator.gene_duplication_rate
cdef double beneficial_mutation_rate = Yeast_Simulator.beneficial_mutation_rate
cdef int n_generations = Yeast_Simulator.number_of_generations
# Cython Imports
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

#FTYPE = np.float
#ctypedef np.float_t dtype_t

##############################################################################################
# Each cell divides at differnt speed. 
##############################################################################################
@cython.boundscheck(False)
def asymmetrical_cell_division(genotypes,cell_groups,PM_fitness_Rate,PM_fitness_Lag,
                                PM,n_generations,n_BN,n_individuals,n_mut,point_mutation,ploidy,text,div_type,s_size,fix_cal,PM_Track,NN,M_RT,LT,M_LT,CRT,CLT,N_HAP1,N_HAP2,N_HAP3,CH,N_CH,MID_Hap):
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
                    genotypes[tp][2] = cell_max_age
                    genotypes[tp][5] = genotypes[k][5]
                    
                    # Copy the Haplotype also
                    CH[tp] = CH[k]
                    N_CH[tp] = N_CH[k] 
                    
                    # Update Rate time of the cell
                    CRT[tp] = CRT[k]
                    
                    # Update Lag time of the cell
                    CLT[tp] = CLT[k]
    
                    PM[tp] = PM[k][:]
                    
                    # Introduce mutations
                    if random.random() <= mutation_rate:
                        
                       # Sending mother and daughter cells to introduce mutations. 
                       genotypes[k],genotypes[tp],PM, CRT,CLT, N_HAP1,N_HAP2,N_HAP3,CH,N_CH,n_mut,MID_Hap = Yeast_Simulator_Mutations.introduce_mutation(genotypes[k],genotypes[tp],
                                                              beneficial_mutation_rate,k,tp,point_mutation,n_mut,PM,PM_fitness_Rate,PM_fitness_Lag,ploidy,fix_cal,CRT,CLT,N_HAP1,N_HAP2,N_HAP3,CH,N_CH,MID_Hap)
            
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
            calculate_freq(CH,i,N_CH,PM,MID_Hap,n_individuals,n_mut)
            (genotypes,cell_groups,PM,CRT,CLT,CH,N_CH) = sampling_individuals(genotypes,n_individuals,i,PM,PM_File,n_BN,LT,M_LT,CRT,CLT,CH,N_CH)
            #calculate_freq(CH,i,N_CH,PM,MID_Hap,n_individuals,n_mut)
        else:
            calculate_freq(CH,i,N_CH,PM,MID_Hap,n_individuals,n_mut)
            (genotypes1,cell_groups,PM,CRT,CLT,CH,N_CH) = sampling_individuals(genotypes,s_size,i,PM,PM_File,n_BN,LT,M_LT,CRT,CLT,CH,N_CH)
            #calculate_freq(CH,i,N_CH,PM,MID_Hap,n_individuals,n_mut)
            sys.exit()
            
            #calculate_mean_cell_division_time(genotypes1,i,n_individuals)
            if fix_cal == 1:
                sys.exit()    
            print "No. of mitoses:",tp
            PM_File.close()
            
            # Fixation analysis: The conditon is to avoid unnecessary fixation calculation where its not needed.
            # Because, Fixation calculation is one of the heavy and time consuming process.
            
            if fix_cal == 1:
                freq_above_0p5,freq_above_0p3,freq_above_0p1 = fixation(i,n_individuals,n_mut,PM_fitness,text)
                
                calculate_haplotype_freq(genotypes1,PM,n_individuals,freq_above_0p5)
            
            return genotypes1,PM,PM_fitness,PM_strand,n_mut

        # Its a small mess up. Have to find a better way.
        groups = []
        groups = cell_groups.keys()
        time = heapq.heappop(groups)

cdef calculate_freq(CH,N_BN,N_CH,PM,MID_Hap,N,n_mut):
    #N = N* 2**5
    #print len(PM)
    #sys.exit()
    freq = {}
    N_WT = 0
    Hap_Freq = {}
    Hap_Freq_N = {}
    N = N * (2**5)
    for i in xrange(N):
        if len(PM[i]) != 0:
            H = "";H_N = 0
            for j in PM[i]:
                if freq.has_key(j):
                    freq[j] += 1
                else:
                    freq[j] = 1
                    
            if len(PM[i]) > 1:
                for j in PM[i]:
                    H +=  MID_Hap[j] + "-" + str(j) + "-"
                if Hap_Freq.has_key(H):
                    Hap_Freq[H] += 1
                else:
                    Hap_Freq[H] = 1
                                 
        else:
            N_WT += 1
  
    #print Hap_Freq
    
    print "Freq:\t",N_BN+1,"\t",0,"\t","WT","\t",N_WT/N,"\t",n_mut
    for i in freq:
        try:
            print "Freq:\t",N_BN+1,"\t",i,"\t",MID_Hap[i],"\t",freq[i]/N,"\t",n_mut
        except:
            pass
    
    print "HFreq:\t",N_BN+1,"\t",0,"\t","WT","\t",N_WT/N
    for i in Hap_Freq:
        try:
            print "HFreq:\t",N_BN+1,"\t",i,"\t",Hap_Freq[i]/N
        except:
            pass
            
    return 0
    
#################################################################################################
# Random Sampling the desired number of individuals after certain number of generations.
#################################################################################################
cdef sampling_individuals(genotypes,int n_individuals,int i,PM1,PM_File,n_BN,LT,M_LT,CRT,CLT,CH,N_CH):

    CRT_Temp = {}
    CLT_Temp = {}
    CH_Temp = {}
    N_CH_Temp = {}
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
    
    """
    NN = 0
    for ii in xrange(len(genotypes)):
        if genotypes[ii][5] == 2:
            NN += 1
            
    print i,"\t",(size-NN)/size,"\t",1
    print i,"\t",NN/size,"\t",2
    """
    
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
        CH_Temp[j] = CH[k]
        N_CH_Temp[j] = N_CH[k]
        
        # Copying PM and GD info of cells
        PM[j] = PM1[k][:]
        
        
        # Write it in a file
        if i == n_BN -1:
            print_in_file(PM_File,PM[j],i)
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
    
    print "CDT\t",NBN+1,"\t",Mean_CLT/len(CLT),"\t","Lag"
    print "CDT\t",NBN+1,"\t",Mean_CRT/len(CRT),"\t","Rate"
    #sys.exit()
    
    return sampled_group,cell_groups,PM,CRT_Temp,CLT_Temp,CH,N_CH


#################################################################################################
# Haplotype Frequency and the CDT for that Haplotype
#################################################################################################
cdef calculate_haplotype_freq(genotypes,PM,n_indi,freq_above_0p5):
    hap_freq = {}; CDT = {};
    # Same Haplotype: different CDT
    SH_D_CDT = 0
    for i in xrange(n_indi):
        key_text = "" 
        for j in xrange(len(PM[i])):
            if j != len(PM[i])-1:
                key_text += str(PM[i][j]) + "-"
            else:
                key_text += str(PM[i][j])
        CDT_temp = genotypes[i][1]
        if hap_freq.has_key(key_text):
            hap_freq[key_text] += 1
            if CDT_temp != CDT[key_text]:
                SH_D_CDT += 1
        else:
            hap_freq[key_text] = 1
            CDT[key_text] = CDT_temp
    print "Number of different Haplotypes:",len(hap_freq)
    print "Same Haplotype with different CDT:",SH_D_CDT
    # Number of distinct Haplotypes based on Mutations above 5% after ASE
    # Distinct Haplotypes
    DH = {}
    for i in hap_freq:
        for j in freq_above_0p5:
            if re.search(str(j),i):
                DH[i] = i
    
    print "Number of mutations above frequency 0p5:",len(freq_above_0p5)                
    print "Number of different haplotypes based on above 0p5:",len(DH)            
        
    k = 0
    for i in hap_freq:
        k += 1
        print "Haplotype\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i]
        
    # Priting the Haplotypes which are having high freq
    for i in hap_freq:
        freq = hap_freq[i]/n_indi
        print "ALL_Haplotype_ASE\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i],"\t",i
        if freq >= 0.01:
            print "S_Haplotype_ASE\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i],"\t",i
        #else:
        #    print "ALL_Haplotype_ASE\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i],"\t",i
            
    return 0
        

#################################################################################################
# Haplotype Frequency and the CDT for that Haplotype
#################################################################################################
def calculate_haplotype_freq_1(genotypes,PM,n_indi,N_BC):
    N_BC += 1
    hap_freq = {}; CDT = {};
    # Same Haplotype: different CDT
    SH_D_CDT = 0
    for i in xrange(n_indi):
        key_text = "" 
        for j in xrange(len(PM[i])):
            if j != len(PM[i])-1:
                key_text += str(PM[i][j]) + "-"
            else:
                key_text += str(PM[i][j])
        CDT_temp = genotypes[i][1]
        if hap_freq.has_key(key_text):
            hap_freq[key_text] += 1
            if CDT_temp != CDT[key_text]:
                SH_D_CDT += 1
        else:
            hap_freq[key_text] = 1
            CDT[key_text] = CDT_temp
            
    #print "Number of different Haplotypes:",len(hap_freq)
    #print "Same Haplotype with different CDT:",SH_D_CDT
        
    k = 0
    for i in hap_freq:
        k += 1
        print "Haplotype\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i]
        
    # Priting the Haplotypes which are having high freq
    text = "S_Haplotype_ABC" + str(N_BC)
    text1 = "ALL_Haplotype_ABC" + str(N_BC)
    for i in hap_freq:
        freq = hap_freq[i]/n_indi
        
        print text1,"\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i],"\t",i
        
        if freq >= 0.01:
            print text,"\t",k,"\t",hap_freq[i]/n_indi,"\t",CDT[i],"\t",i
            
    return 0


cdef print_in_file(FH,PM,i):
    if len(PM) > 0 :
        FH.write(str(i))
        FH.write(",")
        for i in PM:
            FH.write(str(i))
            FH.write(",")
        FH.write("\n")
        
    return 0
    
###############################################################################################    
# This is to calculate the mean cell divison time of the cells at each bottlneck.
###############################################################################################
cdef calculate_mean_cell_division_time(genotypes,int i,int tp):
    # "i" is the bottleneck number.
    cdef int l,n_cells = 0
    cdef float total_time = 0
    cdef float mean_division_time
    
    for l in xrange(tp):
        # Ignore the dead cells.
        if genotypes[l][2] != 0:
            total_time += genotypes[l][1]
            n_cells += 1
        
    mean_division_time = total_time/n_cells
    print i,mean_division_time
    
    return mean_division_time

###############################################################################################
# To find mean number of mutations, gene deletions, gene duplications.
###############################################################################################
cdef calculate_mean_number_of_genetic_variations(genotypes,int tp):
    cdef int i,mean_mutation = 0
    cdef int mean_gene_duplication = 0

    for i in xrange(tp):
        mean_mutation += genotypes[i][3]
        mean_gene_duplication += genotypes[i][4]
    
    print "Mean Mutation:",mean_mutation,tp,(mean_mutation)/float(tp)
    print "Mean Gene Duplication :",mean_gene_duplication,tp,(mean_gene_duplication)/float(tp)
    
    return

##############################################################################################
# To decide whether the cell can divide based on AGE and DIVISION TYPE(either Hap or Dip)
# Its mainly to manage the EXCLUSIVE DIPLOID and EXCLUSIVE HAPLOID division during Backcrossing
##############################################################################################
cdef decide_on_division(genotypes,div_type):
    if genotypes[2] > 0:
        if div_type == 1:
            return 1
        else:
           if genotypes[5] == 2:
              return 1
           else:
              return 0
    else:
        return 0
            
###############################################################################################
# To calculate the Fixaed number of mutations and gene duplications.
# Here we mainly split,extract and format the data
# By calling the cal_fix we calculate the rate of fixation
###############################################################################################
cdef fixation(BN,n,n_mut,PM_fitness,text):
    BN += 1
    
    PM = read_data() 
    
    # Grep Point Mutations and Gene Duplications from last Bottleneck
    PM_fix = greping(BN,PM)

    Last_BN_PM = splitting(PM_fix)
    
    
    U_PM = unique_grep(Last_BN_PM)
    
    freq_above_0p5,freq_above_0p3,freq_above_0p1 = cal_fix(Last_BN_PM,U_PM,n,n_mut,"PM",PM_fitness,text)


    #time_to_fix(PM,BN,n_mut)
    
    return freq_above_0p5,freq_above_0p3,freq_above_0p1

cdef time_to_fix(data,BN,n_mut):
     # Part I: Seperate mutations based on Bottlenecks.
     BN_Mut = []
     U_BN_Mut = []
     for i in xrange(1,BN+1):
        # Grep the data for particular Bottleneck
        temp1 = greping(i,data)
        
        # Splitting the data
        temp2 = splitting(temp1)
        
        # Whole mutation set
        BN_Mut.append(temp2)

        # Unique mutations
        temp3 = unique_grep(temp2)
        U_BN_Mut.append(temp3)
      
     
     """   
     # Part II: Check Each mutation for introduction and extinction time
     for i in range(n_mut):
         t2 = 0; t1 = 0
         for j in range(BN):
             t2 = 0
             for k in range(len(U_BN_Mut[j])):
                 if U_BN_Mut[j][k] == i and t1 == 0:
                    start = j
                    t1 = 1
                 if U_BN_Mut[j][k] == i:
                    t2 = 1
                    
             if t1 == 1 and t2 == 0:
                end = j; time_extinct = (end - start)*5
                print i,"\t",start,"\t",end,"\t",time_extinct
                t2 = 1
                break 
             
         if t1 == 1 and t2 == 1 and j == BN-1:
            end = j;time_extinct = (end - start)*5
            print i,"\t",start,"\t",end,"\t",time_extinct
         
         if t1 == 0 and t2 == 0:
            print  i,"\t0\t0\t0"
     """   
               
###############################################################################################
# To read mutation and gene duplication data
###############################################################################################
cdef read_data():    
    PM1 = open("PM.txt","r")
    PM = PM1.readlines()
    PM1.close()

    return PM
###############################################################################################
# To calculate the Fixated number of mutations and gene duplications.
###############################################################################################
cdef cal_fix(PM,U_PM,n,n_mut,string,fitness,text):
    freq_above_0p5 = {}
    freq_above_0p3 = {}
    freq_above_0p1 = {}
    n_fix = 0 
    for i in U_PM:
         m = 0 
         for j in PM:
            if i == j:
               m += 1
         fix = float(m)/float(n)
         
         if fitness[i] > 0:
            fit_text = "Positive"
         elif fitness[i] == 0:
            fit_text = "Neutral"
         else:
            fit_text = "Deleterious"
        
         if fix >=0.01:
            freq_above_0p1[i] = fix
            
         if fix >=0.03:
            freq_above_0p3[i] = fix
            
         if fix >=0.05:
            freq_above_0p5[i] = fix
            
         print text,"\t",string,"\t",i,"\t",fit_text,"\t",fix,"\t",fitness[i]         
        #print string,"\t",id,"\t",fix,"\t",fitness[id]
         m = 0
         if fix == 1:
            n_fix += 1
    txt = "No. of fixated " + str(string) + ":"        
    print txt,n_fix
    if n_fix != 0:
        print "Rate of fixation:",float(n_fix)/float(n_mut)
        
    return freq_above_0p5,freq_above_0p3,freq_above_0p1   

###############################################################################################
# To extract the unique mutations from the last bottleneck.
###############################################################################################    
cdef unique_grep(PM):
    set = {} 
    map(set.__setitem__, PM, []) 
    
    return set.keys()

###############################################################################################
# To split the mutation numbers and store them in a list from the last bottleneck
###############################################################################################    
cdef splitting(PM):
    Last_BN_PM = []
    for j in PM:       
        temp1 = re.split(",",j)
        #temp2 = re.split("\[",temp1[1])
        #temp3 = re.split("\]",temp2[1])
        #temp4 = re.split("\,",temp3[0])

        for k in xrange(1,len(temp1)-1):    
            Last_BN_PM.append(int(temp1[k]))
        
    return Last_BN_PM
    
###############################################################################################
# To grep the mutation profiles from the last bottleneck
###############################################################################################    
cdef greping(BN,PM):
    BN -= 1
    Last_BN = []
    for i in xrange(len(PM)):
        pattern = "^" + str(BN) + ","
        if re.search(pattern,PM[i]):
           Last_BN.append(PM[i].strip())
    
    return Last_BN
           
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

###############################################################################################
# To group the cels
###############################################################################################
import heapq
def group_the_cells(cell_groups,groups,P_cdt,C_cdt,k,tp):
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
    
    return cell_groups,groups       
          
