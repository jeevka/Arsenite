from __future__ import division
from os import listdir
from os.path import isfile, join
import sys
import re


def store_data(PID,BN,t,R_Time):
    if R_Time.has_key(PID):
        if R_Time[PID].has_key(BN):
            R_Time[PID][BN] += int(float(t))
        else:
            R_Time[PID][BN] = int(float(t))
    else:
        R_Time[PID] = {BN: int(float(t))}
    
    return R_Time

Files = [ f for f in listdir("/home/jeevka/Arsenic_Paper_5/2014_Model_2/100K_Results/") ]

L_Time = {}
R_Time = {}

for i in Files:
    Fname = "100K_Results/" + i
    F = open(Fname,"r")
    for j in F:
        if re.search("Param_ID_",j):
            PID = j.split("_")[2]
            PID = PID.strip()
            
        if re.search("Rate",j): 
            BN = int(j.split()[0])
            t = j.split()[1]
            R_Time = store_data(PID,BN,t,R_Time)
            
        if re.search("Lag",j):
            BN = j.split()[0]
            t = j.split()[1]
            L_Time = store_data(PID,BN,t,L_Time)
    F.close()
    #print R_Time
    #sys.exit()
#print len(R_Time)
#sys.exit()
for i in R_Time:
    for j in R_Time[i]:
        print i,"\t",j,"\t",R_Time[i][j]/5
    #print L_Time