#!/usr/bin/python
import pandas as pd
import numpy as np
import sys
import re
from ast import literal_eval


def returnFitnessPredictCompensationRatio(single_data,feature_name,wt_feature):
    bunch_553_site = literal_eval(single_data["SiteIn553Data_subset"].values[0])
    this_553_site = single_data["SiteIn553Data"].values[0]
    all_n = len(bunch_553_site)
    predict_n = 0
    fitness_diff = literal_eval(single_data["SiteIn553Data"].values[0])
    this_feature_diff = np.abs(AAindex_data.loc[this_553_site,feature_name] - wt_feature)
    bunch_feature_diff  =[np.abs(i-wt_feature) for i in AAindex_data.loc[bunch_553_site,feature_name]]
    
    for i,j  in zip(bunch_diff_feature, fitness_diff):
        if i*j>0:
            predict_n += 1
    return "%d:%d"%(predict_n,all_n)
        
            
featurePredictFit_dir = sys.argv[1]
rescuability_dir = sys.argv[2]
original_dir = sys.argv[3]
aaindex_dir = sys.argv[4]
featurePredictCompen_dir = sys.argv[5]
save_dir = sys.argv[6]


#----------------------------make the rescuability data for subset
featurePredictFit_data = pd.read_csv(featurePredictFit_dir + "FeaturePredictFitness.txt",sep = "\t")
mutation_site82 = featurePredictFit_data["mutation_site"].values
segment_82 = featurePredictFit_data["segment"].values
pick_table_negative_data = pd.read_csv(featurePredictCompen_dir + "pick_table_negtiave.txt",sep = "\t")
Site553_list = [ ]
segment_82_keep = [ ]
for single_mutSite,segment_82_single in zip(mutation_site82,segment_82):
    try:
        pick_table_negagive_data_single = pick_table_negative_data[pick_table_negative_data["mutation_site"] == single_mutSite]["SiteIn553Data"].values[0]
        Site553_list.append(pick_table_negagive_data_single)
        segment_82_keep.append(segment_82_single)
    except IndexError:
        continue

#screen the rescuability data 
#search the site553 in the original data 
for segN in (1,2,3,4,5,6,7,8,10,11,12):
    globals()["AAindex_seg%d"%segN] = pd.read_csv(featurePredictFit_dir + "genotype_553_features_value_segment%d.txt"%segN,sep = "\t")
for segN in (1,2,3,4,5,6,7,8,10,11,12):
    globals()["S_seg%d"%segN] = pd.read_csv(original_dir + "S%d_scaled_info_v2.csv"%segN,sep  = "\t")
for segN in (1,2,3,4,5,6,7,8,10,11,12):
    globals()["Res_seg%d"%segN] = pd.read_csv(rescuability_dir + "rescuability_seg%d.tsv"%segN, sep = "\t")
    
SiteS_list = [ ]
aa_seq_list = [ ] #the aa_seq for the whole genotypes in the file of fitness predict fitness
for i,j in zip(Site553_list, segment_82):
    i_list = literal_eval(i)
    AAindex_data = globals()["AAindex_seg%d"%j]
    SiteS_single = [ ]
    S_seg = globals()["S_seg%d"%j]
    for single553site in i_list:
        aa_seq = AAindex_data.loc[single553site,:]["aa_seq"]
        SiteS_single.append(S_seg[S_seg["aa_seq"] == aa_seq].index[0])
        aa_seq_list.append(aa_seq)
    SiteS_list.append(SiteS_single)
#find the subset for those SiteS_single genotypes!

Site553_site_subset_list = [ ]
Site553_site_list = [ ]
fitness_diff_list = [ ]
for i ,j ,k in zip(SiteS_list, Site553_list,segment_82 ):
    Res_data = globals()["Res_seg%d"%k]
    i_list = i
    j_list = literal_eval(j)
    Site553_site_list.extend(j)
    for singleS_site in i_list:
        fitness_search_list = [ ]
        aaseqSubsetSingle_list = [ ]
        for single_row in range(Res_data.shape[0]):
            SiteInSData_single = literal_eval(Res_data["SiteInSData"].values[single_row])
            if singleS_site in SiteInSData_single:
                aaseq = Res_data["aa_seq"].values[single_row]
                fitness_difflist = literal_eval(Res_data["true_fitness_diff"].values[single_row])
                fitness_diff = fitness_difflist[SiteInSData_single.index[singleS_site]]
                fitness_search_list.append(fitness_diff)
                aaseqSubsetSingle_list.append(aaseq)
        if len(aaseqSubsetSingle_list) == 0:
            Site553_site_subset_list.append([])
        else:
            AAindex_data = globals()["AAindex_seg%d"%k]
            site553_subset_list = [ ]
            for single_aa in aaseqSubsetSingle_list:
                site553_subset_list.append(AAindex_data[AAindex_data["aa_seq"] == single_aa].index[0])
            Site553_site_subset_list.append(site553_subset_list)
            
            fitness_diff_list.append(fitness_search_list)
#make the rescuability data for subset
rescuabilityData_subset = pd.DataFrame({"SiteIn553Data":Site553_site_list,
                                        "SiteIn553Data_subset":Site553_site_subset_list,
                                        "FitnessDiff_minusSubset":fitness_diff_list})
#save the data 
rescuabilityData_subset.to_csv(save_dir + "subset/rescuabilityFileForSubset.tsv",sep = "\t",index = False)

# ----------------------------------repeat the same logic as the superset

rescuabilityData_subset = pd.read_csv(save_dir + "subset/rescuabilityFileForSubset.tsv",sep = "\t")
pick_table_corr_data = pd.read_csv(save_dir + "pick_table_negtiave.txt",sep = "\t")
wt_553_data = pd.read_csv(aaindex_dir + "wild_553_12.txt",sep = "\t")

count_result_list = [ ]
sum_count_list = [ ]
for single_row in range(pick_table_corr_data.shape[0]):
    count_list = [ ]
    feature_name = pick_table_corr_data["feature"].values[single_row]
    SiteIn553Data_list = literal_eval(pick_table_corr_data["SiteIn553Data"].values[single_row]) #compare the features of superset in the AAindex_seg files
    segment_single = pick_table_corr_data["segment"].values[single_row]
    AAindex_data = globals()["AAindex_seg%d"%segment_single]
    wt_feature = wt_553_data.loc[segment_single-1, feature_name]
    for single_site in SiteIn553Data_list:
        try:
            subset_rescuability_single = rescuabilityData_subset[rescuabilityData_subset["SiteIn553Data"] == single_site]
            count_list.append(returnFitnessPredictCompensationRatio(subset_rescuability_single,feature_name,wt_feature))
        except IndexError:
            count_list.append("0:0")
            
    #add into data 
    count_result_list.append(count_list)
    #cycle and summary the count_results
    all_n = 0
    predict_n = 0 
    for i in count_list:
        a,b = i.split(":")
        a = int(a)
        b = int(b)
        all_n += b
        predict_n += a
    sum_count_list.append("%d:%d"%(predict_n,all_n))
    
pick_table_corr_data["sum_count_result"] = sum_count_list
pick_table_corr_data["count_result"] = count_result_list
#save the data  
pick_table_corr_data.to_csv(save_dir + "subset/featurePredictCompensation_subset.tsv",sep = "\t",index = False)