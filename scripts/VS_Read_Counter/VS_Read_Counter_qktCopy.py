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

run_name = sys.argv[1].split("/")[-1]
report = sys.argv[1] + "/" + str(run_name) + ".AssignmentReport"
#print("report=" + report)
file_name = sys.argv[2]
SE_cov_file = sys.argv[1] + "/bbmap_assembly_SE_covstats.txt"
PE_cov_file = sys.argv[1] + "/bbmap_assembly_PE_covstats.txt"

count = 0
join_count = 0
un_count = 0
contig_count = 0
at_organism = False
at_Contig = False
total = 0
print("Sample: " + str(run_name))
print("Category: " + str(sys.argv[2].split("/")[-2]))
from itertools import islice
with open(report) as fin:
    for line in islice(fin, 3, 4):
        print("TotalReads_postQC: " + str(line.split("\t")[1].strip()))

print("Organism\tTotal Read Count\tSingle Read Count\tJoined Read Count\tContig Read Count")
for line in open(file_name):
    if at_organism == True:
        if line.startswith(">"):
            if line.startswith(">un"): #singletons count as 1
                un_count = un_count + 1
                count = count + 1
            elif line.startswith(">join"): #joined reads count as 2
                join_count = join_count + 2
                count = count + 2
            elif line.startswith(">NODE"): #contig:
                if os.system("grep -q "+line.split(">")[1].strip()+" "+PE_cov_file) == True:
                    # print("grep "+line.split(">")[1].strip()+" "+PE_cov_file) # for debugging
                    PE_line = os.system("grep "+line.split(">")[1].strip()+" "+PE_cov_file)
                    PE_total_reads = int(PE_line.split("\t")[6]) + int(PE_line.split("\t")[7])
                    contig_count += PE_total_reads
                    count += PE_total_reads

                    SE_line = os.system("grep "+line.split(">")[1].strip()+" "+SE_cov_file)
                    SE_total_reads = int(SE_line.split("\t")[6]) + int(SE_line.split("\t")[7])
                    contig_count += SE_total_reads
                    count += SE_total_reads
        if line.startswith("#####"):
            print(organism_classifications + "\t" + str(count) + "\t" + str(un_count) + "\t" + str(join_count) + "\t" + str(contig_count))
            at_organism = False
            total = total + join_count + un_count + contig_count
            count = 0
            join_count = 0
            un_count = 0
            contig_count = 0
    if "total number of reads:" in line:
        at_organism = True
        organism_classifications = line.rsplit("\ttotal",1)[0]
	#writing out all levels of classification; this avoids confusion between things such as "Viruses;environmental samples" and "Viruses;dsDNA viruses, no RNA stage;Phycodnaviridae;environmental samples"

print("TOTAL READ COUNT:\t" + str(total))
