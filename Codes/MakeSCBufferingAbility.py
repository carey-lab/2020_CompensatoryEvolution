#!/usr/bin/python
import pandas as pd
import numpy as np
from ast import literal_eval
import matplotlib.pyplot as plt
import math
import sys

compensateAllele_dir = sys.argv[1]
save_dir = sys.argv[2]


#define a function to split the number+alphabet into two:number, alphabet
def split_number_alpha(string):
    number = int(re.findall('\d+',string)[0])
    alphabet = re.findall('\w',string)[-1]
    return number,alphabet
    
def pair_sub_backgroup_data(sub1,sub2,segment):
    aa_seq_contain_sub1_list = []
    aa_seq_contain_sub2_list = []
    #give two list to store the genotype with sub1 and sub2
    for mut_Scer in target_data["whatIwant"]:
        mut_Scer_list = mut_Scer.split(":")
        if sub1 in mut_Scer_list:
            aa_seq_contain_sub1_list.append(mut_Scer)
        elif sub2 in mut_Scer_list:
            aa_seq_contain_sub2_list.append(mut_Scer)    
    #this section find the pair substituions out
    #give two list to store mut_scer
    sub1_mut_scer_keep = []
    sub2_mut_scer_keep = []
    
    for mut_scer_sub1 in aa_seq_contain_sub1_list:
        mut_scer_sub1_list = mut_scer_sub1.split(":")
        #replace the results
        mut_scer_sub2_sub1_list = [sub2 if i==sub1 else i for i in mut_scer_sub1_list]
        mut_scer_sub2_sub1 = ":".join(mut_scer_sub2_sub1_list)
        if mut_scer_sub2_sub1 in aa_seq_contain_sub2_list:
            sub1_mut_scer_keep.append(mut_scer_sub1)
            sub2_mut_scer_keep.append(mut_scer_sub2_sub1)
    
    #return a dataframe
    store_data = pd.DataFrame({sub1:sub1_mut_scer_keep,sub2:sub2_mut_scer_keep})
    return store_data   

#give two aa_seq/mut_list_scer if they are only one distance away, and in the same position
def one_hamming(mut1,mut2):
    mut1_set = set(mut1.split(":"))
    mut2_set = set(mut2.split(":"))
    #there are should be two elements 

    mut1_extra = mut1_set - mut2_set
    mut2_extra = mut2_set - mut1_set
    if len(mut1_extra) != 1:
        return 0
    elif len(mut2_extra) != 1:
        return 0
    else:
        for item in mut1_extra:
            item1_n, item1_al = split_number_alpha(item)
        for item in mut2_extra:
            item2_n, item2_al = split_number_alpha(item)
        if item1_n == item2_n:
            return 1
        else:
            return 0
           
#define functions that make the /return the dataframe which have one distance away from each other
def make_one_distance_data(input_data):
    data_row = input_data.shape[0]
    columns_names_list = list(input_data.columns)
    columns_name1 = columns_names_list[0]
    columns_name2 = columns_names_list[1]

    one_distance_map_col2 = {}
    one_distance_map_col1 = {}
    #using the column2 to search the 1 distance away
    for row in range(data_row-1):
        single_mut_scer_list_col2 = input_data[columns_name2].values[row]
        single_mut_scer_list_col1 = input_data[columns_name1].values[row]
            
        one_distance_map_col2[single_mut_scer_list_col2] = [[]]
        one_distance_map_col1[single_mut_scer_list_col1] = [[]]
        
        second_search_n = row + 1
        for search_row in range(second_search_n , data_row):
            single_search_mut_scer_list_col2  = input_data[columns_name2].values[search_row]
            single_search_mut_scer_list_col1 = input_data[columns_name1].values[search_row]
                
            if one_hamming(single_mut_scer_list_col1 ,single_search_mut_scer_list_col1 ) == 1:
                one_distance_map_col2[single_mut_scer_list_col2][0].append(single_search_mut_scer_list_col2)
                one_distance_map_col1[single_mut_scer_list_col1][0].append(single_search_mut_scer_list_col1)
        
    #section II, using the site row to reapper the detial mut_scer_list original:
    #will return two dataframes
    one_distance_map_col1_data = pd.DataFrame(one_distance_map_col1).T
    one_distance_map_col2_data = pd.DataFrame(one_distance_map_col2).T
    return one_distance_map_col1_data, one_distance_map_col2_data 

#find the difference fitness out 
def return_fitness_difference_return_start(mut_1,mut_2):
    mut_1_fitness = target_data[target_data["whatIwant"] == mut_1]["s"].values[0]
    mut_2_fitness = target_data[target_data["whatIwant"] == mut_2]["s"].values[0]
    return mut_1_fitness, mut_2_fitness, (mut_2_fitness - mut_1_fitness)

#my idea, it is important that to pick all the starting fitness out , or just FE==0 out, anyway, store it 
# this section, try to make a function to pass two sub and return the two lists, and plot
def return_two_sub_one_hamming_fitness_diff_return_start_fitness(sub1,sub2,segment):
    sub12_pair = pair_sub_backgroup_data(sub1,sub2,segment)
    sub1_data, sub2_data = make_one_distance_data(sub12_pair)
    #give two lists to store the fitness difference
    sub1_fitness_diff_list = []
    sub2_fitness_diff_list = []
    start1fitness_list  = []
    end1fitness_list  = []
    start2fitness_list  = []
    end2fitness_list  = []
    
    
    sub1_data_row = sub1_data.shape[0]
    for row in range(sub1_data_row):
        if row%100 == 0:
            print("finish %d"%row)
        sub1_index = sub1_data.index[row]
        sub2_index = sub2_data.index[row]
        for itemsub1 in sub1_data.iloc[row,:].values[0]:
            #print(sub1_index,itemsub1)
            start_fitness, end_fitness, diff_fitness_itemsub1_index = return_fitness_difference_return_start_v1(sub1_index,itemsub1)
            sub1_fitness_diff_list.append(diff_fitness_itemsub1_index)
            start1fitness_list.append(start_fitness)
            end1fitness_list.append(end_fitness)
            
            
            
            
        for itemsub2 in sub2_data.iloc[row,:].values[0]:
            start_fitness, end_fitness,diff_fitness_itemsub2_index = return_fitness_difference_return_start_v1(sub2_index,itemsub2)     
            sub2_fitness_diff_list.append(diff_fitness_itemsub2_index)
            start2fitness_list.append(start_fitness)
            end2fitness_list.append(end_fitness)            
    #return the results
    return sub1_fitness_diff_list, sub2_fitness_diff_list, start1fitness_list,end1fitness_list, start2fitness_list , end2fitness_list
       
#give a function to calculate the skewness of the distribution!
def skewness(distribution):
    s = pd.Series(distribution)
    skew = s.skew()
    return skew


#load the data
Scer_aa_seq_drop_ = "MTEQKALVKRITNETKIQIAISLKGGPLAIEHSIFPEKEAEAVAEQATQSQVINVHTGIGFLDHMIHALAKHSGWSLIVECIGDLHIDDHHTTEDCGIALGQAFKEALGAVRGVKRFGSGFAPLDEALSRAVVDLSNRPYAVVELGLQREKVGDLSCEMIPHFLESFAEASRITLHVDCLRGKNDHHRSESAFKALAVAIREATSPNGTNDVPSTKGVLM"

#-------------------------#

#load the data find pairs of genotypes only with one amino acid difference
sc_data = pd.read_csv(compensateAllele_dir + "SuperCompensation.tsv", sep = "\t")
aaState_list = [ ]
segment_list = [ ]
#find the whole a.a states 
for single_pair,single_seg in zip(sc_data["pair"].values, sc_data["segment"].values):
    if single_seg == 10:
        continue
    stateA,stateB = literal_eval(single_pair)
    if stateA not in aaState_list:
        aaState_list.append(stateA)
        segment_list.append(single_seg)
    if stateB not in aaState_list:
        aaState_list.append(stateB)
        segment_list.append(single_seg)

        
#load the whole 10 original data 
for segN in (1,2,3,4,5,6,7,8,11,12):
    globals()["target_data_seg%d"%segN] = pd.read_csv(compensateAllele_dir + "S%d_scaled_info_v2.csv"%segment_N, sep="\t")

skew_list = []
for aaState,seg in zip(aaState_list,segment_list):
    number_sub,number_al = split_number_alpha(aaState)
    #find the wt amino acids out
    wt_aa = Scer_aa_seq_drop_[number_sub-1]
    target_data = globals()["target_data_seg%d"%seg]
    sub_diff_return_s,wt_diff_return_s,sub_start_s,sub_end_s,wt_start_s,wt_end_s= return_two_sub_one_hamming_fitness_diff_return_start_fitness(aaState,'%d%s'%(number_sub,wt_aa),seg)
    #make it a data
    data_return_start_s = pd.DataFrame({"sub":sub_diff_return_s,"wt":wt_diff_return_s,"sub_start_s":sub_start_s,
                                        "sub_end_s":sub_end_s,"wt_start_s":wt_start_s,"wt_end_s":wt_end_s})
    
    #filter the data
    drop_row_list = [ ]
    for row in range(data_return_start_s.shape[0]):
        sub_start_s = data_return_start_s["sub_start_s"][row]
        wt_start_s = data_return_start_s["wt_start_s"][row]
        sub_end_s = data_return_start_s["sub_end_s"][row]
        wt_end_s = data_return_start_s["wt_end_s"][row]   
        if (sub_start_s ==0) | (wt_start_s ==0) | (sub_end_s == 0) | (wt_end_s == 0):
            drop_row_list.append(row)
    #clean
    data_return_start_s.drop(index=drop_row_list, inplace = True)
    skewness_single = skewness(np.abs(data_return_start_s.iloc[:,0]) - np.abs(data_return_start_s.iloc[:,1]))[0]
    skew_list.append(skewness_single)
#save as a data 
skewness_subsing_data = pd.DataFrame({"sub_name":aaState_list,
                                     "skewness":skew_list})
skewness_subsing_data = skewness_subsing_data[skewness_subsing_data["skewness"] == skewness_subsing_data["skewness"]]
#analysis the compensation ability results
CA_list = [ ]
for singleAA in skewness_subsing_data["sub_name"].values:
    ca_single = 0 
    for single_pair ,single_score in zip(sc_data["pair"].values,sc_data["big05_ratio"].values):
        stateA, stateB = literal_eval(single_pair)
        if singleAA == stateB:
            ca_single += single_score
    CA_list.append(ca_single)    
skewness_subsing_data["compensationAbility"] = CA_list
skewness_subsing_data.to_csv(save_dir + "SubBufferingAbility.tsv",sep = "\t", index = False)
