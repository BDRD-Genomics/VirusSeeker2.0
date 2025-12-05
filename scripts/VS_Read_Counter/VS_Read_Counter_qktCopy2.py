###############################################################################
# Author: Kyle Long <kyle.a.long8.ctr@mail.mil>
# Additional Contact: Regina Cer <regina.z.cer.civ@mail.mil>
#
# License:
# VirusSeeker 2.0 is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3 of the License or any later version. For 
# details, please refer to https://www.gnu.org/licenses/
###############################################################################

import sys
import os
import subprocess as sp

run_name = sys.argv[1].split("/")[-1]
report = sys.argv[1] + "/" + str(run_name) + ".AssignmentReport"
file_name = sys.argv[2]
outfile = open(sys.argv[3],"w")
SE_cov_file = sys.argv[1] + "/bbmap_assembly_SE_covstats.txt"
PE_cov_file = sys.argv[1] + "/bbmap_assembly_PE_covstats.txt"
LR_cov_file = sys.argv[1] + "/minimap_assembly_LR_covstats.txt"

count = 0
join_count = 0
un_count = 0
lr_count = 0
contig_count = 0
at_organism = False
at_Contig = False
total = 0
outfile.write("Sample: " + str(run_name) +"\n")
outfile.write("Category: " + str(sys.argv[2].split("/")[-2]) +"\n")
from itertools import islice
with open(report) as fin:
    for line in islice(fin, 3, 4):
        outfile.write("TotalReads_postQC: " + str(line.split("\t")[1].strip()) +"\n")

outfile.write("Organism\tTotal Read Count\tSingle Read Count\tJoined Read Count\tLong Read Count\tContig Read Count" +"\n")
for line in open(file_name):
    if at_organism == True:
        if line.startswith(">"):
            dashCount = line.count("-")
            noflag = line[1:]
            if line.startswith(">un"): #singletons count as 1
                un_count = un_count + 1
                count = count + 1
            elif line.startswith(">join"): #joined reads count as 2
                join_count = join_count + 2
                count = count + 2
            elif dashCount==4:  #elif line has 4 dashes: #ONT read count as 1
                lr_count += 1
                count += 1

            elif line.startswith(">NODE"): # spades contig:
                cmd = "grep " + line.split(">")[1].strip() + " " + PE_cov_file
                try:
                    PE_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    PE_total_reads = int(PE_line.split("\t")[6]) + int(PE_line.split("\t")[7])
                    contig_count += PE_total_reads
                    count += PE_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue
                cmd = "grep " + line.split(">")[1].strip() + " " + SE_cov_file
                try:
                    SE_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    SE_total_reads = int(SE_line.split("\t")[6]) + int(SE_line.split("\t")[7])
                    contig_count += SE_total_reads
                    count += SE_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue
            elif line.startswith(">contig"): #dragonflye contig
                cmd = "grep " + line.split(">")[1].strip() + " " + PE_cov_file
                try:
                    PE_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    PE_total_reads = int(PE_line.split("\t")[6]) + int(PE_line.split("\t")[7])
                    contig_count += PE_total_reads
                    count += PE_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue
                cmd = "grep " + line.split(">")[1].strip() + " " + SE_cov_file
                try:
                    SE_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    SE_total_reads = int(SE_line.split("\t")[6]) + int(SE_line.split("\t")[7])
                    contig_count += SE_total_reads
                    count += SE_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue
                cmd = "grep " + line.split(">")[1].strip() + " " + LR_cov_file
                try:
                    LR_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    LR_total_reads = int(SE_line.split("\t")[2])
                    contig_count += LR_total_reads
                    count += LR_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue
            elif noflag.isdigit(): #unicycler contig
                cmd = "grep " + line.split(">")[1].strip() + " " + PE_cov_file
                try:
                    PE_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    PE_total_reads = int(PE_line.split("\t")[6]) + int(PE_line.split("\t")[7])
                    contig_count += PE_total_reads
                    count += PE_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue
                cmd = "grep " + line.split(">")[1].strip() + " " + SE_cov_file
                try:
                    SE_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    SE_total_reads = int(SE_line.split("\t")[6]) + int(SE_line.split("\t")[7])
                    contig_count += SE_total_reads
                    count += SE_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue
                cmd = "grep " + line.split(">")[1].strip() + " " + LR_cov_file
                try:
                    LR_line = sp.check_output(cmd, shell=True, encoding="utf-8")
                    LR_total_reads = int(SE_line.split("\t")[2])
                    contig_count += LR_total_reads
                    count += LR_total_reads
                except sp.CalledProcessError: # grep non-zero exit status indicates no match
                    continue

        if line.startswith("#####"):
            outfile.write(organism_classifications + "\t" + str(count) + "\t" + str(un_count) + "\t" + str(join_count) + "\t" + str(lr_count) + "\t" + str(contig_count) + "\n")
            at_organism = False
            total = total + join_count + un_count + contig_count + lr_count
            count = 0
            join_count = 0
            un_count = 0
            contig_count = 0
            lr_count = 0
    if "total number of reads:" in line:
        at_organism = True
        organism_classifications = line.rsplit("\ttotal",1)[0]
	#writing out all levels of classification; this avoids confusion between things such as "Viruses;environmental samples" and "Viruses;dsDNA viruses, no RNA stage;Phycodnaviridae;environmental samples"

outfile.write("TOTAL READ COUNT:\t" + str(total) + "\n")
