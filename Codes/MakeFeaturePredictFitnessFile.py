#!/usr/bin/python
import pandas as pd
import numpy as np
import sys
from scipy import stats
import re

aaindex_dir = sys.argv[1] #the data dir of 553 AAindex features
original_dir = sys.argv[2] #the data dir of S{1-12}_scaled_info_v2.csv except segment9
save_dir = sys.argv[3]

#-----------------------------------the first step is add the 553 aaindex values to the original S data 

aaindex_data = pd.read_csv(aafeature553_dir + "553_features.txt", sep = "\t")
wild_data = pd.read_csv(aafeature553_dir + "wild_553_12.txt",sep = "\t")
def returnFeatureSum(aa_seq,feature_name):
    feature_sum = 0 
    for single_aa in aa_seq:
        feature_sum += aaindex_data[aaindex_data["feature"] == feature_name][single_aa].values[0]
    return feature_sum
        
def calculate_pvalues(df_input):
    df_input = df_input.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df_input.columns)
    for r in df_input.columns:
        dfcols[r] = [round(spearmanr(df_input[r], df_input['fitness'])[1], 4)]
    return dfcols


def calculate_corr(df_input):
    df_input = df_input.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df_input.columns)
    for r in df_input.columns:
        dfcols[r] = [round(spearmanr(df_input[r], df_input['fitness_x'])[0], 4)]
    return dfcols


for segN in (1,2,3,4,5,6,7,8,10,11,12):
    #load and clean the data 
    Original_data = pd.read_csv(original_data_dir+"S%d_scaled_info_v2.csv"%segN,sep="\t")
    #filter the data with sense & size of synonymous substitutions >= 3
    Data_aa_seq_len = [len(i) for i in Original_data.aa_seq]
    Data_aa_seq_nogap = ["_" not in i for i in Original_data.aa_seq ]
    Data = Original_data.loc[(Original_data.nonsense == 0) & (Data_aa_seq_len == stats.mode(Data_aa_seq_len)[0][0]) & Data_aa_seq_nogap]
    Data = Data.loc[(Data["middle"] == 1) & (Data["size"] >= 3)] #remove the WT genotype    
    #cycle the data and add the 553 values
    feature553_list = aaindex_data["feature"].values
    Data_keep = Data.loc[:,["aa_seq","mut_list_Scer","s"]]
    Data_keep.columns = ["aa_seq","mutation_site","fitness"]
    for single_feature in feature553_list:
        feature_sum_list = [ ]
        for single_aaseq in Data["aa_seq"].values:
            feature_sum_list.append(returnFeatureSum(single_aaseq, single_feature))
        Data_keep[single_feature] = feature_sum_list
    #store the data 
    #add the mutation site list
    mutation_site_list = [ ]
    for single_row in range(Data_keep.shape[0]):
        mutation_aa = Data_keep["mutation_site"].values
        mutation_site = [int(i) for i in re.findall("\d+",mutation_aa)]
        mutation_site_list.append(mutation_site)
    Data_keep["mutation_list"] = mutation_site_list
    Data_keep.to_csv(save_dir + "genotype_553_features_value_segment%d.txt"%segN,sep = "\t",index = False)
    
#-------------------------------------the second step is keep the more than 82# unique substituted positions combined

    mutation_clean_list = []
    All_data_row = Data.shape[0]
    for row in range(All_data_row):
        #discard the wild type aa seq
        if Data["dist_Scer"].values[row] == 0:
            print("wild type")
            mutation_clean_list.append([])
            continue     
        else:
            mutation_clean_list.append([int(i) for i in re.findall('(\d+)',All_data["mut_list_Scer"].values[row])])
            
    Data["mutation_list"] = mutation_clean_list
    unique_list = np.unique(mutation_clean_list)
    count_time = []
    for unique_item in unique_list:
        count_time.append(mutation_clean_list.count(unique_item))
    #keep the only genotypes that more than 100
    keep_mutation_site = []
    for i,j in zip(unique_list,count_time):
        if j >= 82:
            keep_mutation_site.append(i)
    #keep the rows that in the keep list
    fitness_list = []
    aa_seq_list = []
    mutation_keep_list = []
    for row in range(All_data_row):
        if Data["mutation_list"].values[row] in keep_mutation_site:
            fitness_list.append(All_data["s"].values[row])
            aa_seq_list.append(All_data["aa_seq"].values[row])
            mutation_keep_list.append(All_data["mutation_list"].values[row])
    All_keep_data = pd.DataFrame({"fitness":fitness_list,"aa_seq":aa_seq_list,"mutation_list":mutation_keep_list})
    #store the All_keep_data
    #merge the aa features data as well as allkeepdata
    AA_features_genotype_data = pd.merge(Data_keep,All_keep_data,on=["aa_seq"])
    #store the merge data
    AA_features_genotype_data.to_csv(save_dir + "553features#82_genotype_seg%d.txt"%segN,sep="\t",index=False)
    
#the next step is make the corraltion of abs(featureOfGenotypes- featureOfWT) and fitness 



segment_list = []
mutation_store_list = []
store_dict = {}
for i in  aaindex_data["feature"].values:
    store_dict[i] = []
store_dict["fitness"] = [ ]

store_dict2= {}
for i in  aaindex_data["feature"].values:
    store_dict2[i] = []
store_dict2["fitness"] = [ ]

for segN in (1,2,3,4,5,6,7,8,10,11,12):
    #load data
    wild_data = pd.read_csv(save_dir + "genotype_553_features_value_segment%d.txt"%segN,sep="\t")
    #extract the wild type things:
    wild_dict = wild_data.iloc[segN-1,:]
    wild_dict["fitness_x"] = wild_type_fitness_dict[segN][1]
    #load aa features
    AA_feature_data = pd.read_csv(save_dir + "553features#82_genotype_seg%d.txt"%segN,sep="\t")
    AA_feature_data_mutation_unique = np.unique(AA_feature_data["mutation_list"].values)
    #keep the unique dataframe
    for unique_mutation_sites in AA_feature_data_mutation_unique:
        segment_list.append(segN)
        mutation_store_list.append(unique_mutation_sites)
        unique_df = AA_feature_data.loc[AA_feature_data["mutation_list"].values == unique_mutation_sites,:]
        unique_df_keep = unique_df.loc[:,list(aaindex_data["feature"].values)]
        unique_df_keep["fitness"] = unique_df["fitness"]
        #combine a new dataframe that can be analysis - wild 
        unique_df3 = pd.DataFrame({})
        for name in unique_df_cal.columns:

            if name != "fitness":
                unique_df3[name] = [np.abs(i) for i in [unique_df_cal[name].values - wild_dict[name]]][0]
            else:
                unique_df3[name] = unique_df_cal[name].values    
        #this step is very important!!,calculate the correlation of dataframe 
        unique_df_corr = calculate_corr(unique_df3)
        for key in store_dict.keys():
            store_dict[key].append(unique_df_corr[key][0])
   
        unique_df_p = calculate_pvalues(unique_df3)
        for key in store_dict2.keys():
            store_dict2[key].append(unique_df_p[key][0])
 
    
#save these two files and filter the correlation to zero if p-value <0.05
store_dict["mutation_site"] = mutation_store_list
store_dict["segment"] = segment_list
df_corr= pd.DataFrame(store_dict)

store_dict2["mutation_site"] = mutation_store_list
store_dict2["segment"] = segment_list
df_p = pd.DataFrame(store_dict2)
keep_column_df = df_p.columns[:-3]
#load 
df_plot2 = df_p.loc[:,keep_column_df]
split_correlation_data = df_corr.iloc[:,0:-3]
select_corr_data = split_correlation_data[df_plot2<0.05]
select_corr_data["mutation_site"] = mutation_store_list
select_corr_data["segment"] = segment_list
#save the data 
select_corr_data.to_csv(save_dir + "FeaturePredictFitness.txt",sep = "\t",index = False)
