#!/usr/bin/python
import pandas as pd
import numpy as np
import sys
import re

rescuabilityFile_dir = sys.argv[1] #the file dir that contain rescuability file
structure_dir = sys.argv[2]  #the file contain the RSA/WCN/Conservation/Interface information for single sites
save_dir = sys.argv[3] #the file to save the rescuability file of each single sites with RSA/WCN/Conservation/Interface information
range_seg_dict = {1:[106,135],2:[136,163],3:[145,170],4:[171,200],5:[181,211],6:[96,125],7:[66,95],8:[56,85],10:[30,47],
                 11:[13,29],12:[6,24]} #the sites of each segment


for segN in (1,2,3,4,5,6,7,8,10,11,12):
    seg_data = pd.read_csv(rescuabilityFile_dir + "rescuability_seg%d.tsv"%segN,sep="\t")
    seg_data = seg_data[(seg_data["s"] == 0) & (seg_data["#genotypesWithAddiSubs"]>=10)] #keep the rescuability for these reliable uunfit genotypes
    globals()["seg%d_data"%segN] = seg_data
    
site_list = []
res_sum = []
count_N = []
res_mean = []
res_sd = []

for site in range(6,212):
    for key in range_seg_dict.keys():
        if (site<=range_seg_dict[key][1]) & (site>=range_seg_dict[key][0]):
            All_data = globals()["seg%d_data"%key]
            All_data_rows = All_data.shape[0]
            ## cycle this segment
            res = 0
            count = 0
            res_single_row_list = []
            for row in range(All_data_rows):
                mutation_site = [int(i) for i in re.findall('(\d+)',All_data["mut_list_Scer"].values[row])]
                if site in mutation_site:
                    res += All_data["rescuability"].values[row]
                    res_single_row_list.append(All_data["rescuability"].values[row])
                    count += 1
            site_list.append(site)
            res_sum.append(res)
            count_N.append(count)
            if count > 0:
                res_mean.append(res/count)
                res_sd.append(np.std(res_single_row_list))
                
            else:
                res_mean.append(0)
                res_sd.append(np.std(res_single_row_list))
    print("end site%d"%site)

singleSiteRes_data = pd.DataFrame({"site":site_list,"#unfitGenotypes":count_N,"rescuability_mean":res_mean,"rescuability_std":res_sd})
#add the information of the WCN-RSA-interface-helix-strand-loop information
singleSiteStru_data = pd.read_csv(structure_dir + "EachSiteInformation.txt",sep = "\t")

singleResAll_data = pd.merge(singleSiteRes_data, singleSiteStru_data)
stru_list = [ ]
for single_row in range(singleResAll_data.shape[0]):
    if singleResAll_data["if_loop"].values[single_row] == 1:
        stru_list.append("loop")
    elif singleResAll_data["if_helix"].values[single_row] == 1:
        stru_list.append("helix")
    else:
        stru_list.append("strand")
singleResAll_data["ss"] = stru_list
singleResAll_data.drop(["if_loop","if_helix","if_strand"],axis = 1,inplace = True)
singleResAll_data.to_csv(save_dir + "singleSiteRescuability.txt",sep = "\t",index = False)
