#!/usr/bin/python
#this script is used to generate the distribution of each amino acid used in each position in His3.
import pandas as pd
import numpy as np
import sys
import os

original_data_dir = sys.argv[1] # the data dir of S{1-12}_scaled_info_v2.csv except segment9
save_dir = sys.arg[2]

site_re_dict = {1:[10,19],2:[9,18],3:[10,18],4:[10,19],5:[10,19],6:[10,19],7:[10,19],8:[10,19],9:[6,11],10:[6,11],11:[6,11],12:[7,12]} #this store 
#the information about mutated positions of each genotype in each segment. For example, 1:[10,19] means that genotypes in segment1 have 
#mutated positions from 0-10 & 19-end.
site_list_dict = [[12,1],[11,1],[12,2],[11,2],
                 [10,1],[9,1],[10,2],[9,2],[8,1],[7,1],[8,2],[7,2],[6,1],[1,1],[6,2],[1,2],[2,1],[3,1],[2,2],[3,2],[4,1],[5,1],[4,2],[5,2]] #this is the connect regions in Hi3. For example, the 1 part in segment12 is adjacent to the 1 part in segment11.

aminoacid_20 = ["G","A","V","L","I","P","F","Y","W","S","T","C","M","N","Q","D","E","K","R","H"]
matrix_6_211 = np.zeros((206,20)) #store the amino acid used frequency in each position of His3.

for segN in range(1,13):
    seg_file = pd.read_csv(original_data_dir + "S%d_scaled_info_v2.csv"%segN,sep="\t")
    #filter the seg_file by milddle and other fators
    Data_aa_seq_len = [len(i) for i in seg_file.aa_seq]
    Data_aa_seq_nogap = ["_" not in i for i in seg_file.aa_seq ]
    seg_file = seg_file.loc[(seg_file.nonsense == 0) & (Data_aa_seq_len == stats.mode(Data_aa_seq_len)[0][0]) & Data_aa_seq_nogap]
    seg_file = seg_file.loc[(seg_file['middle'] == 1)]
    #store the mutated parts of each genotype
    for row in range(seg_file.shape[0]):
        aa_seq  = seg_file['aa_seq'].values[row]
        print(aa_seq[0:site_re_dict[segN][0]]
              , file=open( save_dir + "%d_1.txt"%segN, "a"))
        print(aa_seq[site_re_dict[segN][1]+1:], file=open( save_dir + "%d_2.txt"%segN, "a"))



count = 0
for segN, number in site_list_dict:
    data = pd.read_csv(save_dir + "%d_%d.txt"%(segN,number),names="s",header=None)
    #store the length 
    length = len(data['s'].values[0])
    for f_n in range(length):
        f=lambda x:x[f_n]
        #give a list to store the results
        site_aa = data["s"].apply(f)
        site_aa = list(site_aa)
        
        test_ratio = []
        for aa in aminoacid_20:
            test_ratio.append(list(site_aa).count(aa))
        test_ratio = 100*np.array(test_ratio)/sum(np.array(test_ratio))
        matrix_6_211[count] = test_ratio
        count += 1
    #delete the template data 
    os.remove(save_dir + "%d_%d.txt"%(segN,number))
#save the aa_fre used file for each position
aa_fre = pd.DataFrame(matrix_6_211)
aa_fre.columns = aminoacid_20
aa_fre.to_csv(save_dir +  "aaFrePos6_211.csv",sep=",",index=False)
