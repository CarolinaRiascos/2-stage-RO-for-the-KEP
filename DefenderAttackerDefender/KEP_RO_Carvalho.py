
# Run results: use python os module
# Date: July 7, 2021

# Script to run this executable: nohup python KEP_RO_Carvalho.py 2>&1 &

import os
from pathlib import Path

command = []

Cycle_Length = ["3", "4"]#
Chain_Length = ["3", "4"]#
VertexBudget = ["4"]#, "1", "2",  
ArcBudget = ["0"]#, "1", "2", "3"  
Policy = "Full" 
Formulation = ["Covering", "DoubleCovering"]#  ,"Benders" or "BendersPICEF"
LowerBound = ["Strong"]
WhereToRun = "Server"
TimeLimit = "3600"

InstanceFolder = ["100"]#

exe = "./2StageROThisWork"
LogPrintFolder = "/home/criascos/codes/2StageRO/Output/ROLog"
InstanceName = [[] for i in range(len(InstanceFolder))]

for f in range(len(InstanceFolder)):#Folder
    for s in range(0, 30):#Seeds 
        for n in range(0,10):#NDDs 
            InstanceName[f].append(f'CarvalhoRO2021_{s}_{InstanceFolder[f]}_{n}.txt')
#for f in range(len(InstanceFolder)):
    #print(f'Set {f}: \n')
    #for i in InstanceName[f]:
        #print (i)
    

#for f in Formulation:
for v in VertexBudget:
    for c in Formulation:
        for o in LowerBound:
            for a in ArcBudget:
                for k in Cycle_Length:
                    for l in Chain_Length:
                        for n in range(len(InstanceFolder)):
                            for i in InstanceName[n]:
                                command = command + [ exe + " " + InstanceFolder[n] + " " + i + " " + k + " " + l + " " + v + " " + a + " " + Policy + " " + c + " " + o + " " + WhereToRun + " " + TimeLimit + " > " + LogPrintFolder + "_" + v + "_" + a + "_" + c + "_" + i + "_" + k + "_" + l + ".txt"]
                                #"./2StageROThisWork" "Kidney_Matching_17_1" "KP_Num11_N17_A1.txt" "3" "3" "1" "0" "Full" "Literature" "Server"

for r in range(len(command)):
    print(command[r])
    #os.system(command[r])

