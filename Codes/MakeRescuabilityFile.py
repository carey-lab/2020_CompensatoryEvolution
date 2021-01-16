#!/usr/bin/python
#this script is used to get the rescuability files, 
import pandas as pd
import numpy as np
import sys
from scipy import stats
import re

original_data_dir = sys.argv[1] # the data dir of S{1-12}_scaled_info_v2.csv except segment9
save_dir = sys.argv[2] # the output files of rescuability data for segment from 1-12 except segment9

def BeginBigEnd(single_all_site, original_data, compensated_fitness, aa_seq_compensated):
    compensated_fitness_input = compensated_fitness
    compensated_std = original_data[original_data["aa_seq"] == aa_seq_compensated]["s_std"].values[0]
    
    compensating_fitness = original_data.iloc[single_all_site,:]["s"]
    compensating_std = original_data.iloc[single_all_site,:]["s_std"]
    
    #calculate the region of mean +/- 1.96S.D.
    compensated_fitness_range = [compensated_fitness_input - 1.96*compensated_std,
                                compensated_fitness_input + 1.96*compensated_std]
    compensating_fitness_range = [compensating_fitness - 1.96*compensating_std,
                                 compensating_fitness + 1.96*compensating_std]
  
    if compensating_fitness - compensated_fitness_input >= 0:   
        if compensating_fitness_range[0] >= compensated_fitness_range[1]:  #compensating - 1.96S.D >= compensated + 1.96S.D.
            sig_begin_big_end = compensating_fitness - compensated_fitness_input
        else:
            sig_begin_big_end = 0  
    else:
        sig_begin_big_end = 0
        
    true_fitness_diff = compensating_fitness - compensated_fitness_input
    return sig_true_begin_big_end, true_begin_big_end
    
#calculte the hamming distance of two genotypes
def HammingDistance(char,data_aa,row_Data):
    diff_list = []
    len_char = len(char)
    for i in range(row_Data):
        diff = 0.0
        for j in range(len_char):
            if char[j] != data_aa[i][j]:
                diff += 1.0
        diff_list.append(diff)
    return diff_list

#find the site of Data of substitutions with additional substitutions
def CalculateRescueDist(Unfit_j,Data_dist_j,All_site_editj,edit_dist):
    row_unfit_data = Unfit_j.shape[0]
    row_Data_editj = Data_dist_j.shape[0]
    for i in range(row_unfit_data):
        char = Unfit_j.aa_seq.values[i]
        diff_distance_list = HammingDistance(char,Data_dist_j.aa_seq.values,row_Data_editj)      
        diff_distance_site = [i == edit_dist for i in diff_distance_list]
        All_site_editj.append(list(Data_dist_j.index[diff_distance_site]))


for segN in (1,2,3,4,5,6,7,8,10,11,12):
    
    #---------------------find the sites of genotypes with additional substitutions in Data------------------------
    Original_data = pd.read_csv(original_data_dir+"S%d_scaled_info_v2.csv"%segN,sep="\t")
    #filter the data with sense & size of synonymous substitutions >= 3
    Data_aa_seq_len = [len(i) for i in Original_data.aa_seq]
    Data_aa_seq_nogap = ["_" not in i for i in Original_data.aa_seq ]
    Data = Original_data.loc[(Original_data.nonsense == 0) & (Data_aa_seq_len == stats.mode(Data_aa_seq_len)[0][0]) & Data_aa_seq_nogap]
    Data = Data.loc[(Data["middle"] == 1) & (Data["size"] >= 3) & (Data["dist_Scer"] >= 1)] #remove the WT genotype
    Data = Data.sort_values(by = "dist_Scer")
    Unfit_distance_unique = Data.dist_Scer.sort_values().unique()
    ## here I give 4*25 (25 distance from 1-25) list to store the results of the substitutions sites in Data 
    for edit in range(1,26):
        globals()["All_site_edit%d"%edit] = []
        
    for j in Unfit_distance_unique:
        ## here need know the dist_Scer of each unfit genotype 
        Unfit_j = Data.loc[Data.dist_Scer == j]
        Unfit_j_row = Unfit_j.shape[0]
        ## here fill the column to zero if no such edit distance exist!
        list_zeros = [0 for _ in range(Unfit_j_row)]
        for edit in range(1,26):
            ## the edit must in the Data range(max dist_Scer)
            if (j+edit) <= Data.dist_Scer.max():
                globals()["Data_dist_j%d"%edit] = Data.loc[Data.dist_Scer == (j+edit)]
                CalculateRescueDist(Unfit_j,globals()["Data_dist_j%d"%edit],globals()["All_site_edit%d"%edit],edit)
            else:
                globals()["All_site_edit%d"%edit].extend(list_zeros)

    for edit in range(1,26):
        Data["All_site_%d"%edit] = globals()["All_site_edit%d"%edit]
        
    #-------------------------------make the rescuability data------------------------------------- 
    rescue_data = Data
    Data_row = Data.shape[0]
    #give some lists to store the results
    All_site_list = [ ]
    All_aaSeq_list = [ ]
    All_sSeq_list = [ ]
    All_distScer_list = [ ]
    All_mutListScer_list = [ ]
    All_mutSite_list = [ ]
    Sig_fitnessDiff_list = [ ]
    True_fitnessDiff_list = [ ]
    Rescuability_list = [ ]
    All_N_list = [ ]
    
    for row in range(Data_row):
        all_site_list = [ ]
        sig_fitnessDiff_list = [ ]
        true_fitnessDiff_list = [ ]
        all_n = 0
        for edit in range(1,26):
            if(rescue_data["All_site_edit%d"%edit].values[row] == 0)&(len(rescue_data["All_site_edit%d"%edit].values[row]) > 0):
                continue
            else:
                all_site_list.extend(rescue_data["All_site_edit%d"%edit].values[row])
                all_n +=len(rescue_data["All_site_edit%d"%edit].values[row]) 
       
        if len(all_site_list) > 0:
            mut_site = [int(i) for i in re.findall("\d+",rescue_data["mut_list_Scer"].values[row])]
            All_site_list.append(all_site_list)
            All_N_list.append(all_n)
            All_aaSeq_list.append(rescue_data["aa_seq"].values[row])
            All_sSeq_list.append(rescue_data["s"].values[row])
            All_distScer_list.append(rescue_data["dist_Scer"].values[row])
            All_mutListScer_list.append(rescue_data["mut_list_Scer"].values[row])
            All_mutSite_list.append(mut_site)
            #add the information of fitness difference
            for single_site in all_site_list:
                sig_true_begin_big_end, true_begin_big_end = BeginBigEnd(single_site, original_data, rescue_data["s"].values[row], rescue_data["aa_seq"].values[row])
                sig_fitnessDiff_list.append(sig_true_begin_big_end)
                true_fitnessDiff_list.append(true_begin_big_end)
            
            Sig_fitnessDiff_list.append(sig_fitnessDiff_list)
            True_fitnessDiff_list.append(true_fitnessDiff_list)
        
       #------------------------------------calculate the rescuability-----------------
    
            fitness_big0_sum = 0
            for singleSigFitDiff in Sig_fitnessDiff_list:
                if singleSigFitDiff > 0:
                    fitness_big0_sum += singleSigFitDiff
            Rescuability_list.append(fitness_big0_sum/all_n)
      #---------------------------------------store the data 
     
    Res_data = pd.DataFrame({"aa_seq":All_aaSeq_list, "s":All_sSeq_list, "dist_Scer":All_distScer_list, "mut_list_Scer":All_mutListScer_list,
                             "mut_site":All_mutSite_list, "sig_fitness_diff":Sig_fitnessDiff_list, "true_fitness_diff":True_fitnessDiff_list,
                             "SiteInSData":All_site_list, "#genotypesWithAddiSubs":All_N_list, "Rescuability":Rescuability_list})
    #save the data 
    Res_data.to_csv(save_dir + "rescuability_seg%d.tsv"%segN,index = False,sep = "\t")
    print("End segment%d"%segN)