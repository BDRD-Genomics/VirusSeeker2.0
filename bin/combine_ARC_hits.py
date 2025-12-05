#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:31:04 2024

@author: qthomas
"""

import pandas as pd
import os
import sys, getopt

argv=sys.argv[1:]
opts,args = getopt.getopt(argv,"hr:p:",["rc_path=", "project="])	
#print(opts)
#print(list(sum(opts,())))

project = "out"

for opt,arg in opts:
    if opt == '-h':
        print("combine_ARC_hits.py -r <read_count path> -p <project name>")
        sys.exit()
    elif opt in ("-r", "--rc_path"):
        rc_path = arg
    elif opt in ("-p", "--project"):
        project = arg


#rc_path = "/Users/qthomas/Downloads/VS_processing_hits/ARIA_240611_M07012_VS_outputs/240611_M07012_0141_ARIA_Hyb_VS_counts/VS_processed/"
#project = "ARIA_240611_M07012"
#/Users/qthomas/Downloads/AFI_GEO_202405/VS_results/processed_ARCS/AFI_GEO_USAMRDG_KMED-053-1_2022_no_hits.txt

pattern="_AccurateReadCounts_normalized.txt"


arcs_tmp = os.popen("ls "+rc_path+"/*"+pattern).read().split('\n')
arc_files = [x for x in arcs_tmp if x]

#print(pos_arc_files[0])
#print(os.path.getsize(pos_arc_files[0]))
#tmp = pd.read_csv(pos_arc_files[0], sep='\t')

#pos and neg files: Family, Beyond_Family, total_counts, normalized_counts
# no hit files: Family, Beyond_Family, total_counts, normalized_counts, Full_lineage

main_frame = pd.DataFrame(columns=["Full_Lineage","Family", "Beyond_Family"])
#main_frame['total_counts']=main_frame['total_counts'].astype(float)
#main_frame['normalized_counts']=main_frame['normalized_counts'].astype(float)
#main_frame.total_counts.astype(float)
#main_frame.normalized_counts.astype(float)
#main_frame = pd.DataFrame()
main_frame_total = main_frame
main_frame_normalized = main_frame

empty_files=[]
all_samples = []

for file in arc_files:
    #print(file)
    sample =  os.path.basename(file).replace("_AccurateReadCounts_normalized.txt","")
    all_samples+=[sample]
    empty_check = os.path.getsize(file)
    #print(sample)
    if (empty_check == 0):
        empty_files+=[sample]
    else:
        tmp_df = pd.read_csv(file, sep='\t', skiprows=4, names=["Full_Lineage","Total_Read_Count","Single_Read_Count","Joined_Read_Count","Long_Read_Count","Contig_Read_Count","Family", "Beyond_Family", "normalized_family", "normalized_rpm"])
        print(tmp_df.head())
        tmp_df.drop(tmp_df[tmp_df.Full_Lineage == "TOTAL READ COUNT:"].index, inplace=True)
        tmp_df.replace("No family found", None, regex=True, inplace=True)
        tmp_df.Beyond_Family.fillna('', inplace=True)
        tmp_df.Beyond_Family.astype(str)
        num_rows = len(tmp_df)
        #print(num_rows)
        if num_rows != 0:
            #print(tmp_df["Beyond_Family"])
            tmp_df[sample+"_total"] = tmp_df["Total_Read_Count"].astype(float)
            tmp_df[sample+"_normalizedFamily"] = tmp_df["normalized_family"].astype(float)
            tmp_df[sample+"_normalizedRPM"] = tmp_df["normalized_rpm"].astype(float)
            subdf_normalizedFam = tmp_df[["Full_Lineage","Family", "Beyond_Family", sample+"_normalizedFamily"]]
            subdf_normalizedRPM = tmp_df[["Full_Lineage","Family", "Beyond_Family", sample+"_normalizedRPM"]]
            subdf_total = tmp_df[["Full_Lineage","Family", "Beyond_Family",sample+"_total"]]

            #print()
            #print(tmp_df.head())
            normalized_main_frame = pd.merge(main_frame_normalized, subdf_normalizedFam, how="outer", on=["Full_Lineage", "Family", "Beyond_Family"])
            normalized_main_frame_rpm = pd.merge(main_frame_normalized, subdf_normalizedRPM, how="outer", on=["Full_Lineage", "Family", "Beyond_Family"])
            total_main_frame = pd.merge(main_frame_total, subdf_total, how="outer", on=["Full_Lineage","Family", "Beyond_Family"])
            main_frame_total = total_main_frame
            main_frame_normalized = normalized_main_frame
            main_frame_normalized_rpm = normalized_main_frame_rpm
        else:
            empty_files+=[sample]


df_save_no_hits_total = main_frame_total
df_save_no_hits_normFam = main_frame_normalized
df_save_no_hits_normRPM = main_frame_normalized_rpm

curr_cols = list(df_save_no_hits_total.columns)
#print(curr_cols)

print("Total number of samples: ", len(all_samples))


#print(empty_files)
add_cols = list(set([x for x in empty_files if empty_files.count(x)==3]))
#print(add_cols)
for sam in add_cols:
    main_frame_total[sam] = None
    main_frame_normalized[sam] = None

#print(main_frame_total.head())
main_frame_total.to_csv(rc_path+project+"_total_hits_merged.csv")
main_frame_normalized.to_csv(rc_path+project+"_normalized_family_hits_merged.csv")
main_frame_normalized_rpm.to_csv(rc_path+project+"_normalized_RPM_hits_merged.csv")







