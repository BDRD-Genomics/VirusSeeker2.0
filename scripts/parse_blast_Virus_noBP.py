#python parse blast using pandas

# 1) read blast output table
# 2) remove duplicate rows
# 3) groupby query ID and target ID
# 4) combine multiple query/target results into a single line
    # This means calculate new alignment length, new coverage coordinates and count the number of times this query aligned to this target
# 5) Add lineage information
# 6) Sort df by query and then by bitscore

###########################################################################
# Data Prep  #
###########################################################################


import sys
import os
import pandas as pd
import numpy as np
import argparse
import subprocess as sp
import sqlite3
from ete3 import NCBITaxa
import configparser
from Bio import SeqIO

# Get paths from config file
config = configparser.ConfigParser()
config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'VS.cfg')
config.read(config_path)
paths = config['paths']
vhunter = paths['vhunter']
ncbidb = paths['ncbi_taxadb']
ncbi = NCBITaxa(dbfile = ncbidb)

# open a connection to database
print("Connecting to vhunter db...")
try:
	connector = sqlite3.connect(vhunter)
except:
	print("Cannot connect to DB: "+vhunter)

cursor_dbhsqlite = connector.cursor()
print("Connected")

parser = argparse.ArgumentParser(
                prog='Parse MMseqs Outputs',
                description='',
                epilog='Parse the blast table from sequence comparative algorithm (MMSeqs,blastn,blastx)',
                usage="parse_blast_noBP.py -i <input dir> -f <input file> -o <output dir>")
# add options
parser.add_argument("-i", "--input") # blast/mmseqs input table
parser.add_argument("-f", "--fasta") # blast/mmseqs input table
parser.add_argument("-d", "--dir") # blast/mmseqs input table

args = parser.parse_args()
print(args.input)
if (args.input == None):
        print(parser.usage)
        exit(0)
bdir = args.dir if args.dir else "."
if not os.path.exists(bdir):
    os.mkdir(bdir)
fasta_input = args.fasta
blastout=args.input

# get file name prefix
file_prefix_arr = blastout.split(".")
print(file_prefix_arr)
file_name_prefix = ".".join(file_prefix_arr[0:3])
#file_name_suffix = "."+str(file_prefix_arr[-1])
if file_prefix_arr[3] == "mmseqsv":
	outFile = bdir+"/"+file_name_prefix+".mmseqsv.parsed"
	e_cutoff = 1e-5
	extra_rem = False
	out_suffix = ".MMseqs_VIRUSDB"
	os.system('sed -i "$ d" {0}'.format(bdir+"/"+blastout))
	vhunter_nuc = True
	BX=False
elif file_prefix_arr[3] == "BLASTX_VIRUSDB":
	outFile = bdir+"/"+file_name_prefix+".BLASTX_VIRUSDB.parsed"
	e_cutoff = 1e-3
	extra_rem = True
	out_suffix = ".BXVIRUSDB"
	vhunter_nuc = False
	BX=True

#try:
#    out = open(outFile, 'w')
#except IOError:
#	print("can not open file "+outFile)

phage_list = ["Uroviricota", "Loebvirae","Trapavirae","Sangervirae","Caudovirales","Caudoviricetes","unclassified bacterial viruses"]
phage_lowercase=[x.lower() for x in phage_list]

unassigned = [] # query should be kept for further analysis
virusseq = []
phageseq = []

total_records = 0

print("parsing blast output files...\n\n")

input_file = bdir+"/"+blastout
#custom_fields=["qseqid","sacc","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sseqid","qlen","slen"]


###########################################################################
# Functions  #
###########################################################################

def read_FASTA_data(fastaFile):
	fa_dict = SeqIO.index(fastaFile, "fasta")
	return fa_dict

def get_lineage(hit_name, description):
    #print(hit_name)
    #print(description)
    #hit_name=row["sseqid"]
    #description=row["sallgi"]
    if BX:
        sth = cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_prot where accession_version = '"+hit_name+"'")
    else:
        sth = cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_nucl where accession_version = '"+hit_name+"'")
    ref = sth.fetchone()
    if ref: # some gi don't have record in gi_taxid_nucl
        taxID = ref[2]
        taxon_name = ncbi.get_taxid_translator([taxID])
        formatted_lineage=""
        if taxon_name:
            lineage = ncbi.get_lineage(taxID)
            lineage.pop(0)
            if lineage:
                for temp_node_id in lineage:
                     temp_name = ncbi.get_taxid_translator([temp_node_id])
                     formatted_lineage += temp_name[temp_node_id]+";"
        else:
             formatted_lineage="undefined taxon: "+description
    else:
        formatted_lineage="undefined taxon: "+description
    return formatted_lineage

def getAssignment(formatted_lineage):
    assignment =""
    if "Metazoa" and "Homo" in formatted_lineage and formatted_lineage.startswith("cellular_organisms;Eukaryota"):
        assignment="Human"
    elif "Metazoa" and "Mus" in formatted_lineage and formatted_lineage.startswith("cellular_organisms;Eukaryota"):
        assignment="Mouse"
    elif "Metazoa" in formatted_lineage and formatted_lineage.startswith("cellular_organisms;Eukaryota"):
        assignment="Other"
    elif formatted_lineage.startswith("cellular_organisms;Bacteria"):
        assignment="Bacteria"
    elif formatted_lineage.startswith("cellular_organisms;Eukaryota;Opisthokonta;Fungi"):
        assignment="Fungi"
    elif formatted_lineage.startswith("Viruses"):
        assignment= "Phage" if any(phage in formatted_lineage.lower() for phage in phage_lowercase) else "Virus"
    elif formatted_lineage.startswith("undefined taxon"):
        assignment="unassigned"
    else:
        assignment="Other"
    return assignment

def get_final_df(testdf,vseq,pseq,unassi):
    concatenated_df = pd.DataFrame()
    for unique_id,unique_df in testdf.groupby("qseqid", as_index=False):
        #print(unique_id)
        #print(unique_df)
        row_to_return = pd.Series([])
        num_assignments=len(unique_df['test'].iloc[0])
        if num_assignments == 1:
            if unique_df['test'].iloc[0][0] == "Virus":
                vseq+=[unique_id]
                #print(unique_df['bitscore'].idxmax())
                row_to_return = unique_df.loc[unique_df['bitscore'].idxmax()]
                #print(vseq)
            elif unique_df['test'].iloc[0][0] == "Phage":
                pseq+=[unique_id]
                row_to_return = unique_df.loc[unique_df['bitscore'].idxmax()]
                #print(pseq)
            else:
                unassi+=[unique_id]
                #print(unassi)
        elif num_assignments > 1:
            if "Virus" in unique_df['test'].iloc[0]:
                if num_assignments ==2 and "unassigned" in unique_df['test'].iloc[0]:
                    row_to_return = tmpdf.loc[tmpdf['bitscore'].idxmax()]
                    vseq+=[unique_id]
                else:
                    tmpdf = unique_df[unique_df["assignment"] != "Virus"]
                    if len(np.unique(tmpdf.test))==1:
                        row_to_return = tmpdf.loc[tmpdf['bitscore'].idxmax()]
                        if row_to_return["assignment"] == "Phage":
                            pseq+=[unique_id]
                        else:
                            unassi+=[unique_id]
                    #print(row_to_return)
                    else:
                        unassi+=[unique_id]
            else:
                unassi+=[unique_id]
        if len(row_to_return) > 0:
            concatenated_df = pd.concat([concatenated_df,row_to_return.to_frame().T])
    print(concatenated_df.head())
    return concatenated_df,virusseq,phageseq,unassigned
   #check to see if length of test col is greater than 1
    #if len is one, find the row with teh highest bitscore 
    #if len is greater than one, check to see which taxa are listed
        #if virus and *OTHER* then remove rows where assignment=virus
        #if *OTHER* and *OTHER* then


###################################################################################

###########################################################################
# Step 1:  #
###########################################################################

custom_fields=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sallgi","qlen","slen"]

# Read blast-like table
blast_df = pd.read_csv(args.input, delimiter="\t", low_memory=True, engine="pyarrow", names=custom_fields)
#print(str(blast_df.shape[0]))
blast_df = blast_df.drop_duplicates()
#print(blast_df.head())


# Select the rows using the obtained indices
#filter the dataframe so that all of the single query-subject alignment are condensed into a single line
blast_df = blast_df[blast_df['evalue'] < e_cutoff]
blast_df["aln_count"] = blast_df.groupby(["qseqid", "sseqid"])["sseqid"].transform('size')
#find the row with the highest bit score for each query-subject pair
idx = blast_df.groupby(["qseqid", "sseqid"])['bitscore'].idxmax()
sub_blast_df1 = blast_df.loc[idx]

#remove rows below the evalue cutoff
#sub_blast_df1 = sub_blast_df[sub_blast_df['evalue'] < e_cutoff]

print("Adding lineage...")
#sub_blast_df1["lineage"] = sub_blast_df1.apply(lambda row: get_lineage(row["sseqid"], row["sallgi"]), axis=1)
sub_blast_df1["lineage"] = sub_blast_df1.apply(lambda row: get_lineage(row["sseqid"], row["sallgi"]), axis=1)
print("Adding assignment...")
sub_blast_df1["assignment"] = sub_blast_df1.apply(lambda row: getAssignment(row["lineage"]), axis=1)
#print(sub_blast_df1.head(25))

#Next step!!!!!!!!!!!!!
#summarize assignments per query
print("Summarizing assignments")
test_col = sub_blast_df1.groupby("qseqid")['assignment'].apply(lambda x: list(set(x))).reset_index(name="test")
#print(sub_blast_df1.head(30))
testdf = pd.merge(sub_blast_df1, test_col, on="qseqid", how="left")
print(testdf.head(25))
#tt= testdf.set_index("qseqid").groupby("qseqid", as_index=False).apply(lambda sub: filter_bad_rows(sub))
final_blast_df,virusseq,phageseq,unassigned  = get_final_df(testdf, virusseq,phageseq,unassigned)

#print(virusseq[0])
## format df output
list_of_col_names = ["qseqid", "qlen", "assignment", "lineage", "sseqid", "slen", "sallgi", "length", "pident_formatted", "target_span1", "target_span2", "evalue"]
#reformat percent identity with %
final_blast_df['pident_formatted'] = final_blast_df['pident'].astype(str) + "%"
#reformat sstart and send for aln span
final_blast_df['target_span1'] = "["+final_blast_df['sstart'].astype(str)
final_blast_df['target_span2'] = final_blast_df['send'].astype(str)+ "]"


formatted_blast_df = final_blast_df[list_of_col_names]
print(formatted_blast_df.head())

'''
18 rows: 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
       'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sallgi', 'qlen',
       'slen', 'aln_count', 'lineage', 'assignment'
'''
#print("Filtering lineage...")
#sub_blast_df1["lin_count"] = sub_blast_df1.groupby(["qseqid", "lineage"])["lineage"].transform('size')
#print(len(sub_blast_df1))
#idx2 = sub_blast_df1.groupby(["qseqid", "lineage"])['bitscore'].idxmax()
#sub_blast_df2 = sub_blast_df1.loc[idx2]
#print(len(sub_blast_df2))
#print(sub_blast_df2.columns)


# print(sub_blast_df1.head(n=15))
testdf.to_csv("test_full.parsed", index=False, sep='\t', header=False)
print("Writing output to: "+outFile)
formatted_blast_df.to_csv(outFile, index=False, sep='\t', header=False)


print("Reading fasta file")
#file = bdir+"/"+file_name_prefix+".fasta"
unassigned2 = []
seq_dict = read_FASTA_data(fasta_input)
seq_headers = list(seq_dict.keys())

hit_headers = set(virusseq + phageseq)
unassigned2 = list(set(seq_headers) ^ hit_headers)

outFile = bdir+"/"+file_name_prefix+out_suffix+"_hit.fa"
try:
	out2 = open(outFile, "w")
except:
	print("Can't create output file : "+outFile)
print("VS list length " + str(len(virusseq)))
for seq_name in virusseq:
	out2.write(">"+seq_name+"\n")
	seq=seq_dict[seq_name].seq
	out2.write(str(seq)+"\n")
out2.close()

outFile = bdir+"/"+file_name_prefix+out_suffix+"_phage.fa"
try:
	out2 = open(outFile, "w")
except:
	print("Can't create output file : "+outFile)

print("Phage list length "+str(len(phageseq)))
for seq_name in phageseq:
	out2.write(">"+seq_name+"\n")
	seq=seq_dict[seq_name].seq
	out2.write(str(seq)+"\n")
out2.close()

outFile = bdir+"/"+file_name_prefix+out_suffix+"_filtered.fa"
try:
	out2 = open(outFile, "w")
except:
	print("Can't create output file : "+outFile)

print("Unassigned list length "+str(len(unassigned2)))
for seq_name in unassigned2:
	out2.write(">"+seq_name+"\n")
	seq=seq_dict[seq_name].seq
	out2.write(str(seq)+"\n")
out2.close()



