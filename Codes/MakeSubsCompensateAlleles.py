#!/usr/bin/python
import pandas as pd
import numpy as np
from ast import literal_eval
import matplotlib.pyplot as plt
import math
import sys

compensate_dir = sys.argv[1] #the rescuability files
original_dir = sys.argv[2] #the S_{1-12} files
save_dir = sys.argv[3]

for segN in (1,2,3,4,5,6,7,8,10,11,12):
    #-------------------- store the genotypes with one additional substitution and fitness 
    
    original_data = pd.read_csv(original_dir + "S%d_scaled_info_v2.csv"%segN,sep="\t")
    fitness_diff_list = []
    all_site_list = []

    compen_data = pd.read_csv(compensate_dir + "rescuability_seg%d.csv"%segN,sep="\t")
    compen_data = compen_data[compen_data["dist_Scer"] != 0] #remove WT genotype
    compen_data_row = compen_data.shape[0]
    for row in range(compen_data_row):
        fitness_diff_single_list = []
        all_site_this_list = []


        site_list = eval(compen_data["SiteInSData"].values[row])
        fitness_list = eval(compen_data["true_fitness_diff"].values[row])
        this_dist_Scer = compen_data["dist_Scer"].values[row]
        for single_original_site, single_fitness_diff in zip(site_list,fitness_list):
            superset_dist_Scer = original_data.iloc[single_original_site,:]["dist_Scer"]
            if (super_dist_Scer - this_dist_Scer) != 1:
                continue
            else:

                fitness_diff_single_list.append(single_fitness_diff)
                all_site_this_list.append(single_original_site)


        fitness_diff_list.append(fitness_diff_single_list)
        all_site_list.append(all_site_this_list)


    compen_data["SiteInSData"] = all_site_all_list
    compen_data["true_fitness_diff"] = fitness_diff_all_list
    #clean the data 
    compen_data_clean = compen_data[compen_data["SiteInSData"] != '[]']
    #----------------------------------find the details of the additional one substitutions
    
    mut_list_Scer_all_list = []
    mut_list_Scer_all_superset_list = []
    data_compen_row = compen_data.shape[0]
    for row in range(data_compen_row):
        superset_add_list = []
        mut_list_Scer_single_list = (compen_data_clean["mut_list_Scer"][row]).split(":")
        all_single_site_list = compen_data_clean["SiteInSData"][row]
        for all_single_site in all_single_site_list:
            mut_list_Scer_superset_list =(original_data.iloc[all_single_site,:]["mut_list_Scer"]).split(":")
            mut_list_Scer_superset_add = list(set(mut_list_Scer_super_list) - set(mut_list_Scer_single_list))[0]
            superset_add_list.append(mut_list_Scer_super_add)
        mut_list_Scer_all_list.append(mut_list_Scer_single_list)
        mut_list_Scer_all_superset_list.append(super_add_list)
    compenData_WithOneAddSub = pd.DataFrame({"oriset":mut_list_Scer_all_list,"addset":mut_list_Scer_all_super_list,"fit_diff":compen_data_clean["fitness_list"]})
    #----------------------------------- make the whole possible combination of the original mutation set and one additional substitution

    compenData_WithOneAddSub_row = compenData_WithOneAddSub.shape[0]
    AddOneCompen_list = []
    for single_row in range(compenData_WithOneAddSub_row):
        oriset_list = compenData_WithOneAddSub["oriset"][single_row]
        addset_list = compenData_WithOneAddSub["addset"][single_row]
        c = itertools.product(oriset_list,addset_list)
        AddOneCompenaSingle_list = list(c)
        AddOneCompen_list.append(AddOneCompenaSingle_list)
    compenData_WithOneAddSub["pair_list"] = AddOneCompen_list
    
    
    #-----------------------------------------find the pair sites for further analysis
    
    pair_all_list = []
    #give a list to store all pairs
    compenData_WithOneAddSub_row = compenData_WithOneAddSub.shape[0]
    for single_row in range(compenData_WithOneAddSub_row):
        pair_single_list = compenData_WithOneAddSub["pair_list"][single_row]
        pair_all_list.extend(pair_single_list)
    #set and reomove the duplication
    unique_pair_list = set(pair_all_list)
    #store the site of each mut pair
    mut_pair_map = {}
    for mut_pair in unique_pair_list:
        mut_pair_map[mut_pair] = [[]]     
    #search the site in this step    
    for mut_pair in mut_pair_map.keys():
        for row in range(compenData_WithOneAddSub_row):
            pair_single_list = eval(compenData_WithOneAddSub_row["pair_list"][row])
            if mut_pair in pair_single_list:
                mut_pair_map[mut_pair][0].append(row)
    #save as a dataframe
    str_mut_pair = {}
    for mut_pair in mut_pair_map:
        str_mut_pair[str(mut_pair)] = mut_pair_map[mut_pair]
    #save as a datafrme
    mut_pair_site_data = pd.DataFrame(str_mut_pair).T
    #--------------------------------------------- make the final data 
    mut_pair_site_data.columns = ["pair_site"]
    mutPair_data_row = mut_pair_site_data.shape[0]
    diff_all_list = []
    for row in range(mutPair_data_row):
        diff_single_list = []
        index_single = mut_pair_site_data.index[row]
        nothing,add_muat = index_single
        
        single_pair_site = mut_pair_site_data["pair_site"][row]
        for single_site in single_pair_site:
            fit_diff_list = compenData_WithOneAddSub.iloc[single_site,:]["fit_diff"]
            mutate_list = compenData_WithOneAddSub.iloc[single_site,:]["addset"]
            fitness_diff = fit_diff_list[mutate_list.index(add_muat)]
            diff_single_list.append(fitness_diff)
        diff_all_list.append(diff_single_list)
    mut_pair_site_data["diff_fitness"] = diff_all_list
    #---------------------count the %of genotypes which has fitness difference>0.5
    count_N_list = []
    data_row = mut_pair_site_data.shape[0]
    for row in range(data_row):
        first_N = 0
        biger_05_N = 0 
        diff_list = mut_pair_site_data["diff_fitness"][row]
        for diff_fitness in diff_list:
            first_N += 1
            if diff_fitness >= 0.5:
                biger_05_N +=  1
        count_N_list.append((first_N,biger_05_N))
    mut_pair_site_data["all_big05"] = count_N_list
    #save the data 
    mut_pair_site_data.to_csv(save_dir + "CompensatePair_seg%d.tsv"%segN,sep="\t")

#-------------------------merge the 12 segments

for segN in (1,2,3,4,5,6,7,8,10,11,12):  
    mut_pair_site_data  = pd.read_csv(save_dir + "CompensatePair_seg%d.tsv"%segN,sep="\t")
    globals()["compenAllele_%d"%segN] = mut_pair_site_data
    globals()["compenAllele_%d"%segN]["seg"] = [segN]*len(globals()["compenAllele_%d"%segN])

compenAllele_all = pd.concat([compenAllele_1,compenAllele_2,compenAllele_3,compenAllele_4,
                         compenAllele_5,compenAllele_6,compenAllele_7,
                              compenAllele_8,compenAllele_10,compenAllele_11,
                              compenAllele_12],axis=0,ignore_index=True)
#also add the ratio
ratio_list = []
for row in range(compenAllele_all.shape[0]):
    count_tuple = eval(compenAllele_all["all_big05"][row])
    a,b = count_tuple
    ratio_list.append(b/a)
compenAllele_all["bigg05_ratio"] = ratio_list
compenAllele_all.columns = ["pair","add_mutation_site","fitness_diff","all_big05",
                       "segment",
                       "big05_ratio"]
#drop the name 
compenAllele_all.drop(["add_mutation_site"],axis = 1,inplace = True)
compenAllele_all.to_csv(save_dir + "SuperCompensation.tsv", sep="\t",index=False)