import pandas as pd
import numpy as np
import sys
from ast import literal_eval

res_dir = sys.argv[1] #the dir containing the rescuability file
stru_dir = sys.argv[2] #the dir containing the structure information file
aafeature553_dir = sys.argv[3]
save_dir = sys.argv[4] #the dir to save the data 

#----functions----#
def returnStruInfo(mut_list, feature_name):
    feature_list = [ ]
    mut_list_eval = literal_eval(mut_list)
    for singleMutSite in mut_list_eval:
        feature_list.append(struc_data[struct_data["site"] == singleMutSite][feature_name].values[0])
    return np.max(feature_list), np.min(feature_list), np.mean(feature_list), np.median(feature_list)
def returnSSInfo(mut_list,feature_name):
    feature_list = [ ]
    mut_list_eval = literal_eval(mut_list)
    for singleMutSite in mut_list_eval:
        feature_list.append(struc_data[struct_data["site"] == singleMutSite][feature_name].values[0])    
    feature_ratio = 100*np.sum(feature_list)/len(feature_list)
    return feature_ratio
def returnAAFeatureValue(aa_seq,feature_name,segment):
    feature_list = aaindex_data[aaindex_data["feature"] == "ANDN920101"]
    feature_sum = 0
    absfeatureDiff = 0 
    for i in aa_seq:
        feature_sum += feature_list[i].values[0]
    absfeatureDiff = np.abs(feature_sum - wild_aaindex_data[feature_name][segment-1])
    return feature_sum, absfeatureDiff
    

      
#load datas 
struc_data = pd.read_csv(stru_dir + "EachSiteInformation.txt",sep = "\t")
aaindex_data = pd.read_csv(aafeature553_dir + "553_features.txt", sep = "\t")
wild_aaindex_data = pd.read_csv(aafeature553_dir + "wild_553_12.txt",sep = "\t")
#filter threshold:only make the machine learning data for unfit genotypes with #genotypes with addtional substitutions >= 10
for segN in (1,2,3,4,5,6,7,8,10,11,12):
    seg_data = pd.read_csv(data_dir + "rescuability_seg%d.tsv"%segN,sep="\t")
    seg_data = seg_data[(seg_data["s"] == 0) & ( seg_data["All_N"]>=10)]
    globals()["seg%d_data"%segN] = seg_data
#merge the segdata into the whole data
all_segment_list = []
all_segment_list.extend([1]*seg1_data.shape[0])
all_segment_list.extend([2]*seg2_data.shape[0])
all_segment_list.extend([3]*seg3_data.shape[0])
all_segment_list.extend([4]*seg4_data.shape[0])
all_segment_list.extend([5]*seg5_data.shape[0])

all_segment_list.extend([6]*seg6_data.shape[0])
all_segment_list.extend([7]*seg7_data.shape[0])
all_segment_list.extend([8]*seg8_data.shape[0])
all_segment_list.extend([10]*seg10_data.shape[0])
all_segment_list.extend([11]*seg11_data.shape[0])
all_segment_list.extend([12]*seg12_data.shape[0])
res_data = pd.concat([seg1_data,seg2_data,seg3_data,
                     seg4_data,seg5_data,seg6_data,seg7_data,seg8_data,seg10_data,
                     seg11_data,seg12_data],axis=0,ignore_index=True)
res_data["segment"] = all_segment_list   
#add the information about WCN-RSA-Inteface-helix-strand-loop-conservation-information
res_data["ln_allN"] = np.log10(res_data["#genotypesWithAddiSubs"])
stru_info_list = ["WCN","conservation","RSA"]
ss_info_list = ["loop","helix","strand","interface"]
for singleStruFeature in stru_info_list: #add the information of structure 
    max_list = [ ]
    min_list = [ ]
    mean_list = [ ]
    median_list = [ ]
    sum_list = [ ]
    for single_row in range(res_data.shape[0]):
        i,j,m,n,q = returnStruInfo(res_data["mut_site"].values[single_row], singleStruFeature)  
        max_list.append(i)
        min_list.append(j)
        mean_list.append(m)
        median_list.append(n)
        sum_list.append(q)
        
    res_data["max_%s"%singleStruFeature] = max_list
    res_data["min_%s"%singleStruFeature] = min_list
    res_data["mean_%s"%singleStruFeature] = mean_list
    res_data["median_%s"%singleStruFeature] = median_list
    
for singleSSFeature in ss_info_list:#add the information about secondray structure
    ratio_list = [ ]
    for single_row in range(res_data.shape[0]):
        ratio_list.append(res_data["mut_site"].values[single_row], singleSSFeature)
    res_data["%s_ratio"%singleSSFeature] = ratio_list
                                

feature_name_list = wild_aaindex_data.columns #add the 553 features information
for singleFeature in feature_name_list:
    feature_sum_list = [ ]
    feature_diff_list = [ ]
    for single_aa in zip(res_data["aa_seq"].values, res_data["segment"].values):
        i,j = returnAAFeatureValue(single_aa,singleFeature,segment)
        feature_sum_list.append(i)
        feature_diff_list.append(j)
    res_data["%d_true"%singleFeature] = feature_sum_list
    res_data["%d_diff"%singleFeature] = feature_diff_list

 #save the data
res_data.drop(["s","mut_list_Scer","mut_site","sig_fitness_diff","true_fitness_diff","SiteInSData"],axis = 1,inplace = True)
res_data.to_csv(save_dir + "MLPrepare_data.tsv", sep = "\t",index = False)
            
    
                                
                                
    
