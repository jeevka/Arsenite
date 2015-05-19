from __future__ import division
import sys
import re


F1 = open(sys.argv[1],"r")
F2 = F1.readlines()

F3 = open(sys.argv[2],"r")
F4 = F3.readlines()

Param_ID = int(sys.argv[3])

COF = float(sys.argv[4])

F5 = open("Intro_Freq.csv","a")

Freq = {}
Fitness = {}
for i in F2:
    if re.search("^ASE",i):
        ID = int(re.split("\t",i)[2])
        Fre = re.split("\t",i)[4]
        Fit = re.split("\t",i)[5]
        if float(Fre) >= COF:
            Freq[ID] = Fre
            Fitness[ID] = Fit
    if re.search("No. of fixated",i):
        break


for i in Freq:
    F5.write(str(Param_ID))
    F5.write("\t")
    F5.write(str(Freq[i]))
    F5.write("\t")
    F5.write(str(Fitness[i].strip()))
    F5.write("\t")
    F5.write(str(F4[i-1].strip()))
    F5.write("\n")
F5.close()