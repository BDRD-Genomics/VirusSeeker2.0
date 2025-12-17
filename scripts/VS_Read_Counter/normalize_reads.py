from ete3 import NCBITaxa
import sys, getopt
import os
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE
ncbi = NCBITaxa(dbfile = "/export/database/taxonomy/taxa.sqlite")

argv=sys.argv[1:]
opts,args = getopt.getopt(argv,"hi:",["in="])
for opt,arg in opts:
    if opt == 'h':
        print("normalize_reads.py -i <inputfile>")
        sys.exit()
    elif opt in ("-i", "--in"):
        ARC_file = arg

file_name = ARC_file.rsplit(".",1)[0]
print(file_name)

# This file comes from NCBI's viruses.tsv (ftp:/genomes/GENOME_REPORTS) a few additional family sizes have been added due to incomplete naming and other issues. Up to date Oct2023
virus_genomes_file = "/export/virusseeker/scripts/qkt_scripts/viralFamily_genomeSize.txt"
vf_dict = pd.read_csv(virus_genomes_file, sep="\t", header=0, names=["fam","size"]).set_index('fam')['size'].to_dict()
#print(vf_dict)

output = open(file_name+"_normalized.txt", "w")

with open (ARC_file) as ARC:
    for i, line in enumerate(ARC):
        output.write(line)
        if i==2:
            total_QC_reads = int(line.split(": ")[1])
            break
output.close()

ARC_df = pd.read_csv(ARC_file, delimiter='\t', skiprows=3) # read AccurateReadCounts into df
last_row = ARC_df.tail(1).copy()
ARC_df.drop(ARC_df.tail(1).index, inplace=True) # remove last row which is the total counts

family_arr = []
normalized_reads_arr = []
sub_family_arr = []
rpm_arr=[]

def normalization(counts, QC_reads, genome_length):
    return ((counts/QC_reads)/genome_length)*1000
def rpm_norm(counts, QC_reads):
    return (counts/QC_reads)*1000000

for ind in ARC_df.index:
    full_taxid, total_counts = ARC_df['Organism'][ind], ARC_df['Total Read Count'][ind]
    parsed = full_taxid.split(';')
    #print(full_taxid)
    parsed.pop()
    
    taxid=ncbi.get_name_translator([parsed[-1]])
    lineage=ncbi.get_lineage(taxid.get(parsed[-1])[0])
    rank = ncbi.get_rank(lineage)

    family_key=list(rank.keys())
    family_value=list(rank.values())
    #print(family_key)
    try:
        family_idx=family_value.index("family")
        family_dict=ncbi.get_taxid_translator([family_key[family_idx]])
        family=family_dict.get(family_key[family_idx])
        #print(family)
        family_arr.append(family)
        command_str = "cat "+virus_genomes_file+" | grep "+family
        try:
            taxid_family_index = parsed.index(family)
            #print(full_taxid[taxid_family_index:])
            sub_family = ';'.join(parsed[taxid_family_index+1:])
            sub_family_arr.append(sub_family)
            #print(command_str)
        except ValueError:
            sub_family_arr.append('N/A')

    except ValueError:
        print("This entry has no family order")
        norm_rpm = rpm_norm(total_counts,total_QC_reads)
        rpm_arr.append(norm_rpm)
        family_arr.append("No family found")
        normalized_reads_arr.append("No family found")
        sub_family_arr.append('N/A')
        continue
    try:
        '''
        with Popen(command_str, shell=True, stdout=PIPE) as process:
            family_genomes_df = pd.read_csv(process.stdout, delimiter='\t', header=None)
        #print(type(family_genomes_df))
        #print(type(family_genomes_df[6]))
        if family_genomes_df.dtypes[6] == 'object':
             family_genomes_df[7] = family_genomes_df.apply(lambda r: r[6] if type(r[6])==float else np.nan, axis=1)
             family_genomes_df.dropna(inplace=True)
        #print(family_genomes_df[6].astype(float).mean())
        '''        
        mean_genome_family = vf_dict[family]
        #print(mean_genome_family)
        normalized_counts = normalization(total_counts,total_QC_reads,mean_genome_family)
        #print(normalized_counts)
        normalized_reads_arr.append(normalized_counts)
        norm_rpm = rpm_norm(total_counts,total_QC_reads)
        rpm_arr.append(norm_rpm)
    except pd.errors.EmptyDataError:
        #print("No family genomes found.")
        norm_rpm = rpm_norm(total_counts,total_QC_reads)
        rpm_arr.append(norm_rpm)
        normalized_reads_arr.append("No family genomes found")
        continue

#print(ARC_df['Organism'].tolist())
#print(family_arr)
#for i1, i2 in zip(ARC_df['Organism'].tolist(), family_arr):
#	print(i1, i2)

ARC_df['Family'] = family_arr
ARC_df['below Family'] = sub_family_arr
ARC_df['normalized_counts'] = normalized_reads_arr
ARC_df['normalized_rpm'] = rpm_arr
last_row['normalized_rpm'] = 'NaN'
last_row['Family'] = 'NaN'
last_row['normalized_counts'] = 'NaN'
last_row['below Family'] = 'NaN'
final_df = pd.concat([ARC_df,last_row])


final_df.to_csv(file_name+"_normalized.txt", mode='a', sep='\t', index=False, header=True)

