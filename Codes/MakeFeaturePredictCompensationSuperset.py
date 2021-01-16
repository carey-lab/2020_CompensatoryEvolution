#!/usr/bin/python
import pandas as pd
import numpy as np
import sys
import re
from ast import literal_eval


def returnFitnessPredictCompensationRatio(single_data,feature_name,wt_feature):
    bunch_553_site = literal_eval(single_data["553_site_bunch"].values[0])
    this_553_site = single_data["553_site"].values[0]
    all_n = len(bunch_553_site)
    predict_n = 0
    fitness_diff = literal_eval(single_data["fitness_diff"].values[0])
    this_feature_diff = np.abs(AAindex_data.loc[this_553_site,feature_name] - wt_feature)
    bunch_feature_diff  =[np.abs(i-wt_feature) for i in AAindex_data.loc[bunch_553_site,feature_name]]
    
    for i,j  in zip(bunch_diff_feature, fitness_diff):
        if i*j<0:
            predict_n += 1
    return "%d:%d"%(predict_n,all_n)
        
            
featurePredictFit_dir = sys.argv[1]
rescuability_dir = sys.argv[2]
original_dir = sys.argv[3]
aaindex_dir = sys.argv[4]
save_dir = sys.argv[5]

#------------make the pick table datas 
featurePredictFit_data = pd.read_csv(featurePredictFit_dir + "FeaturePredictFitness.txt",sep = "\t")

corr_list = [ ]
feature_name_list = [ ]
segment_list = [ ]
mutation_site_list = [ ]
for single_row in range(featurePredictFit_data.shape[0]):
    for single_col in range(553):
        feature_name = featurePredictFit_data.columns[single_col]
        corr_single = featurePredictFit_data.iloc[single_row,single_col]
        mutation_site_single = featurePredictFit_data.loc[single_row,"mutation_site"]
        segment_single = featurePredictFit_data.loc[single_row,"segment"]
        if corr_single < 0:
            feature_name_list.append(feature_name)
            segment_list.append(segment_single)
            mutation_site_list.append(mutation_site_single)
            corr_list.append(corr_single)
#make and store as a dataframe
pick_table_corr_negative_data = pd.DataFrame({"correlation":corr_list, "segment":segment_list, "mutation_site":mutation_site_list,
                                              "feature":feature_name_list})

#add the information about sites in 553 AAfeature genotype data 
for segN in (1,2,3,4,5,6,7,8,10,11,12):
    globals()["AAindex_seg%d"%segN] = pd.read_csv(featurePredictFit_dir + "genotype_553_features_value_segment%d.txt"%segN,sep = "\t")
    
siteIn553_data_list  = [ ]
for single_row in range(pick_table_corr_negative_data.shape[0]):
    mutation_site_single  =  pick_table_corr_negative_data["mutation_site"].values[single_row]
    segment_single = pick_table_corr_negative_data["segment"].values[single_row]
    mutation_site_single_data = globals()["AAindex_seg%d"%segment_single][globals()["AAindex_seg%d"%segment_single]["mutation_list"]=mutation_site_single]
    site_index_list = list(mutation_site_single_data.index)
    siteIn553_data_list.append(site_index_list)
pick_table_corr_negative_data["SiteIn553Data"] =  siteIn553_data_list


#------------------make additional files that serve as the feature predict compensation based on rescuabilityt data 
#first load the original data 
for segN in (1,2,3,4,5,6,7,8,10,11,12):
    globals()["S_seg%d"%segN] = pd.read_csv(original_dir + "S%d_scaled_info_v2.txt"%segN,sep = "\t")
    
unique_sites_combination = featurePredictFit_data["mutation_site"].values
aa_seq_list = [ ]
fitness_diff_list = [ ]
single_553_site_list = [ ]
bunch_553_site_list = [ ]
mut_sites_list = [ ]
for segN in (1,2,3,4,5,6,7,8,10,11,12):
    rescuability_data = pd.read_csv(rescuability_dir + "rescuability_seg%d.tsv"%segN, sep = "\t")
    for single_row in range(rescuability_data.shape[0]):
        mut_site_single = rescuability_data["mut_site"].values[single_row]
        if mut_site_single in unique_sites_combination:
            mut_sites_list.append(mut_site_single)
            aa_seq_list.append(rescuability_data["aa_seq"].values[single_row])
            fitness_diff_list.append(rescuability_data["true_fitness_diff"].values[single_row])
            single_553_site = globals()["AAindex_seg%d"%segN][globals()["AAindex_seg%d"%segN]["aa_seq"] == rescuability_data["aa_seq"].values[single_row]].index[0]
            single_553_site_list.append(single_553_site)
            #searth the 553 sites for genotypes with additional substitutions
            bunch_original_site = literal_eval(rescuability_data["SiteInSData"].values[single_row])
            bunch_sites_single_list = [ ]
            for single_site in bunch_original_site:
                aa_seq_single = globals()["S_seg%d"%segN].loc[single_site,:]["aa_seq"]
                site_553_single =globals()["AAindex_seg%d"%segN][globals()["AAindex_seg%d"%segN]["aa_seq"] == rescuability_data["aa_seq"].values[single_row]].index[0]
                bunch_sites_single_list.append(site_553_single)
            bunch_553_site_list.append(bunch_sites_single_list)
#make the serve data 
rescuability_featurePredictCompensation_data = pd.DataFrame({"aa_seq":aa_seq_list, "mut_site":mut_sites_list,
                                                             "553_site":single_553_site_list,"fitnessDiff_additionsubs":fitness_diff_list,
                                                             "553SiteBunch_additionsubs":bunch_553_site_list})

#sava both data 
pick_table_corr_negative_data.to_csv(save_dir + "pick_table_negtiave.txt",sep = "\t",index = False)
rescuability_featurePredictCompensation_data.to_csv(save_dir + "superset/rescuabilityFileForSuperset.txt", sep = "\t",index = False)

#------------find the number of superset genotypes of which feature predict compensation
pick_table_corr_data = pd.read_csv(save_dir + "pick_table_negtiave.txt",sep = "\t")
rescuability_superset_data = pd.read_csv(save_dir + "superset/rescuabilityFileForSuperset.txt", sep = "\t")
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
            superset_rescuability_single = rescuability_superset_data[rescuability_superset_data["553_site"] == single_site]
            count_list.append(returnFitnessPredictCompensationRatio(superset_rescuability_single,feature_name,wt_feature))
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
pick_table_corr_data.to_csv(save_dir + "superset/featurePredictCompensation_superset.tsv",sep = "\t",index = False)