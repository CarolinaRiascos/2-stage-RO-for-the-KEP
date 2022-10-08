
# Run results: use python os module
# Date: July 7, 2021

# Script to run this executable: nohup python KEP_RO.py 2>&1 &

import os
from pathlib import Path

command = []

Cycle_Length = ["3"]#, "4"
Chain_Length = ["3"]#, "4"
VertexBudget = ["3"]#,  
ArcBudget = ["0"]#, "1", "2", "3"  
Policy = "Full" #seconds
Formulation = ["Covering"]# , "Literature"
WhereToRun = "Server"
TimeLimit = "3600"

InstanceFolder = ["Kidney_Matching_17_1", "Kidney_Matching_67_3"]#, "Kidney_Matching_70_6",  "Kidney_Matching_73_9", "Kidney_Matching_134_6", "Kidney_Matching_18_2",  "Kidney_Matching_33_1", "Kidney_Matching_35_3", "Kidney_Matching_36_4", "Kidney_Matching_140_12", "Kidney_Matching_147_19","Kidney_Matching_268_12", "Kidney_Matching_281_25", "Kidney_Matching_294_38", "Kidney_Matching_537_25", "Kidney_Matching_563_51", "Kidney_Matching_588_76", "Kidney_Matching_1075_51", "Kidney_Matching_1126_102","Kidney_Matching_1177_153","Kidney_Matching_2150_102", "Kidney_Matching_2252_204","Kidney_Matching_2355_307"

InstanceName = [
    ["KP_Num11_N17_A1.txt","KP_Num12_N17_A1.txt","KP_Num13_N17_A1.txt","KP_Num14_N17_A1.txt","KP_Num15_N17_A1.txt","KP_Num16_N17_A1.txt","KP_Num17_N17_A1.txt","KP_Num18_N17_A1.txt","KP_Num19_N17_A1.txt","KP_Num20_N17_A1.txt"],
    #["KP_Num21_N18_A2.txt","KP_Num22_N18_A2.txt","KP_Num23_N18_A2.txt","KP_Num24_N18_A2.txt","KP_Num25_N18_A2.txt","KP_Num26_N18_A2.txt","KP_Num27_N18_A2.txt","KP_Num28_N18_A2.txt","KP_Num29_N18_A2.txt","KP_Num30_N18_A2.txt"],
    #["KP_Num41_N33_A1.txt","KP_Num42_N33_A1.txt","KP_Num43_N33_A1.txt","KP_Num44_N33_A1.txt","KP_Num45_N33_A1.txt","KP_Num46_N33_A1.txt","KP_Num47_N33_A1.txt","KP_Num48_N33_A1.txt","KP_Num49_N33_A1.txt","KP_Num50_N33_A1.txt"],
    #["KP_Num51_N35_A3.txt","KP_Num52_N35_A3.txt","KP_Num53_N35_A3.txt","KP_Num54_N35_A3.txt","KP_Num55_N35_A3.txt","KP_Num56_N35_A3.txt","KP_Num57_N35_A3.txt","KP_Num58_N35_A3.txt","KP_Num59_N35_A3.txt","KP_Num60_N35_A3.txt"],
    #["KP_Num61_N36_A4.txt","KP_Num62_N36_A4.txt","KP_Num63_N36_A4.txt","KP_Num64_N36_A4.txt","KP_Num65_N36_A4.txt","KP_Num66_N36_A4.txt","KP_Num67_N36_A4.txt","KP_Num68_N36_A4.txt","KP_Num69_N36_A4.txt","KP_Num70_N36_A4.txt"],
    ["KP_Num81_N67_A3.txt","KP_Num82_N67_A3.txt","KP_Num83_N67_A3.txt","KP_Num84_N67_A3.txt","KP_Num85_N67_A3.txt","KP_Num86_N67_A3.txt","KP_Num87_N67_A3.txt","KP_Num88_N67_A3.txt","KP_Num89_N67_A3.txt","KP_Num90_N67_A3.txt"],
    #["KP_Num100_N70_A6.txt","KP_Num91_N70_A6.txt","KP_Num92_N70_A6.txt","KP_Num93_N70_A6.txt","KP_Num94_N70_A6.txt","KP_Num95_N70_A6.txt","KP_Num96_N70_A6.txt","KP_Num97_N70_A6.txt","KP_Num98_N70_A6.txt","KP_Num99_N70_A6.txt"],
    #["KP_Num101_N73_A9.txt","KP_Num102_N73_A9.txt","KP_Num103_N73_A9.txt","KP_Num104_N73_A9.txt","KP_Num105_N73_A9.txt","KP_Num106_N73_A9.txt","KP_Num107_N73_A9.txt","KP_Num108_N73_A9.txt","KP_Num109_N73_A9.txt","KP_Num110_N73_A9.txt"],
    #["KP_Num121_N134_A6.txt","KP_Num122_N134_A6.txt","KP_Num123_N134_A6.txt","KP_Num124_N134_A6.txt","KP_Num125_N134_A6.txt","KP_Num126_N134_A6.txt","KP_Num127_N134_A6.txt","KP_Num128_N134_A6.txt","KP_Num129_N134_A6.txt","KP_Num130_N134_A6.txt"]
    #["KP_Num131_N140_A12.txt","KP_Num132_N140_A12.txt","KP_Num133_N140_A12.txt","KP_Num134_N140_A12.txt","KP_Num135_N140_A12.txt","KP_Num136_N140_A12.txt","KP_Num137_N140_A12.txt","KP_Num138_N140_A12.txt","KP_Num139_N140_A12.txt","KP_Num140_N140_A12.txt"],
    #["KP_Num141_N147_A19.txt","KP_Num142_N147_A19.txt","KP_Num143_N147_A19.txt","KP_Num144_N147_A19.txt","KP_Num145_N147_A19.txt","KP_Num146_N147_A19.txt","KP_Num147_N147_A19.txt","KP_Num148_N147_A19.txt","KP_Num149_N147_A19.txt","KP_Num150_N147_A19.txt"],#
    #["KP_Num161_N268_A12.txt","KP_Num162_N268_A12.txt","KP_Num163_N268_A12.txt","KP_Num164_N268_A12.txt","KP_Num165_N268_A12.txt","KP_Num166_N268_A12.txt","KP_Num167_N268_A12.txt","KP_Num168_N268_A12.txt","KP_Num169_N268_A12.txt","KP_Num170_N268_A12.txt"],
    #["KP_Num171_N281_A25.txt","KP_Num172_N281_A25.txt","KP_Num173_N281_A25.txt","KP_Num174_N281_A25.txt","KP_Num175_N281_A25.txt","KP_Num176_N281_A25.txt","KP_Num177_N281_A25.txt","KP_Num178_N281_A25.txt","KP_Num179_N281_A25.txt","KP_Num180_N281_A25.txt"],
    #["KP_Num181_N294_A38.txt","KP_Num182_N294_A38.txt","KP_Num183_N294_A38.txt","KP_Num184_N294_A38.txt","KP_Num185_N294_A38.txt", "KP_Num186_N294_A38.txt","KP_Num187_N294_A38.txt","KP_Num188_N294_A38.txt","KP_Num189_N294_A38.txt","KP_Num190_N294_A38.txt"]#
    #["KP_Num201_N537_A25.txt","KP_Num202_N537_A25.txt","KP_Num203_N537_A25.txt","KP_Num204_N537_A25.txt","KP_Num205_N537_A25.txt","KP_Num206_N537_A25.txt","KP_Num207_N537_A25.txt","KP_Num208_N537_A25.txt","KP_Num209_N537_A25.txt","KP_Num210_N537_A25.txt"],
    #["KP_Num211_N563_A51.txt","KP_Num212_N563_A51.txt","KP_Num213_N563_A51.txt","KP_Num214_N563_A51.txt","KP_Num215_N563_A51.txt","KP_Num216_N563_A51.txt","KP_Num217_N563_A51.txt","KP_Num218_N563_A51.txt","KP_Num219_N563_A51.txt","KP_Num220_N563_A51.txt"],
    #["KP_Num221_N588_A76.txt","KP_Num222_N588_A76.txt","KP_Num223_N588_A76.txt","KP_Num224_N588_A76.txt","KP_Num225_N588_A76.txt","KP_Num226_N588_A76.txt","KP_Num227_N588_A76.txt","KP_Num228_N588_A76.txt","KP_Num229_N588_A76.txt","KP_Num230_N588_A76.txt"],#
    #["KP_Num241_N1075_A51.txt","KP_Num242_N1075_A51.txt","KP_Num243_N1075_A51.txt","KP_Num244_N1075_A51.txt","KP_Num245_N1075_A51.txt","KP_Num246_N1075_A51.txt","KP_Num247_N1075_A51.txt","KP_Num248_N1075_A51.txt","KP_Num249_N1075_A51.txt","KP_Num250_N1075_A51.txt"],
    #["KP_Num251_N1126_A102.txt","KP_Num252_N1126_A102.txt","KP_Num253_N1126_A102.txt","KP_Num254_N1126_A102.txt","KP_Num255_N1126_A102.txt","KP_Num256_N1126_A102.txt","KP_Num257_N1126_A102.txt","KP_Num258_N1126_A102.txt","KP_Num259_N1126_A102.txt","KP_Num260_N1126_A102.txt"],
    #["KP_Num261_N1177_A153.txt","KP_Num262_N1177_A153.txt","KP_Num263_N1177_A153.txt","KP_Num266_N1177_A153.txt","KP_Num267_N1177_A153.txt","KP_Num268_N1177_A153.txt","KP_Num269_N1177_A153.txt","KP_Num270_N1177_A153.txt","KP_Num264_N1177_A153.txt","KP_Num265_N1177_A153.txt"],#
    #["KP_Num281_N2150_A102.txt","KP_Num282_N2150_A102.txt","KP_Num283_N2150_A102.txt","KP_Num284_N2150_A102.txt","KP_Num285_N2150_A102.txt","KP_Num286_N2150_A102.txt","KP_Num287_N2150_A102.txt","KP_Num288_N2150_A102.txt","KP_Num289_N2150_A102.txt","KP_Num290_N2150_A102.txt"],
    #["KP_Num291_N2252_A204.txt","KP_Num292_N2252_A204.txt","KP_Num293_N2252_A204.txt","KP_Num294_N2252_A204.txt","KP_Num295_N2252_A204.txt","KP_Num296_N2252_A204.txt","KP_Num297_N2252_A204.txt","KP_Num298_N2252_A204.txt","KP_Num299_N2252_A204.txt","KP_Num300_N2252_A204.txt"],
    #["KP_Num301_N2355_A307.txt","KP_Num302_N2355_A307.txt","KP_Num303_N2355_A307.txt" ,"KP_Num304_N2355_A307.txt", "KP_Num305_N2355_A307.txt", "KP_Num306_N2355_A307.txt","KP_Num307_N2355_A307.txt","KP_Num308_N2355_A307.txt","KP_Num309_N2355_A307.txt","KP_Num310_N2355_A307.txt"]#
   ] #


exe = "./2StageROThisWork"
LogPrintFolder = "/home/criascos/codes/2StageRO/Output/ROLog"
FileNames = []

#for f in Formulation:
for v in VertexBudget:
    for c in Formulation:
        for a in ArcBudget:
            for k in Cycle_Length:
                for l in Chain_Length:
                    for n in range(len(InstanceFolder)):
                        for i in InstanceName[n]:
                            FileNames.append(i)
                            command = command + [ exe + " " + InstanceFolder[n] + " " + i + " " + k + " " + l + " " + v + " " + a + " " + Policy + " " + c + " " + WhereToRun + " " + TimeLimit]#+ " > " + LogPrintFolder + "_" + c + "_" + i + "_K" + k + "_L" + l + ".txt"
                            "./2StageROThisWork" "Kidney_Matching_17_1" "KP_Num11_N17_A1.txt" "3" "3" "1" "0" "Full" "Literature" "Server"

for r in range(len(command)):
    print(command[r])
    os.system(command[r])

