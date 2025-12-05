from Bio import SearchIO
from Bio import SeqIO
from ete3 import NCBITaxa

import sqlite3
import sys
import os
import configparser


Usage = '''
This script accepts a BLASTn output file and parse the information \n
python script <dir><blast output file> \
<dir> = directory that blast output file resides in, without last "/"
<blast output file> = name of the BLASTx output file'''

try:
    bdir = sys.argv[1]
    blastout = sys.argv[2]
except:
	sys.exit(1)

# Get paths from config file
config = configparser.ConfigParser()
config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'VS.cfg')
config.read(config_path)
paths = config['paths']
vhunter = paths['vhunter']
ncbidb = paths['ncbi_taxadb']
ncbi = NCBITaxa(dbfile = ncbidb)

# open a connection to database
try:
	connector = sqlite3.connect(vhunter)
except:
	print("Cannot connect to DB: "+vhunter)

cursor_dbhsqlite = connector.cursor()

def read_FASTA_data(fastaFile):
	fa_dict = SeqIO.index(fastaFile, "fasta")
	return fa_dict

# function to determine the taxonomy lineage for a given blast hit
def PhyloType(lineage_ref, result_ref, hit_ref):
	assignment_ref = {}
	assigned = 0
	description = ""
	lineage = ""
	#This for loop basically just grabs the scientific name for all the taxids in the lineage and saves it to a single variable
	for temp_node_id in lineage_ref:
		temp_name = ncbi.get_taxid_translator([temp_node_id])
		lineage += temp_name[temp_node_id]+";"
		#print(lineage)

	# check to see if it is a human sequence
	#if (scalar @{$lineage_ref} >= 4) {
	if len(lineage_ref) >= 4:
		node_id = lineage_ref[3]
		obj = ncbi.get_taxid_translator([node_id])
		name = obj[node_id]
		#print(name)
		if name == "Metazoa":
			# make assignment
			for temp_node_id in lineage_ref:
				temp_obj = ncbi.get_taxid_translator([temp_node_id])
				temp_name = temp_obj[temp_node_id] # double check this attribute
				if temp_name == "Homo":
					if "Homo" not in assignment_ref.keys():
						target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
						Homo_desc = "\t".join(["Homo",hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
						assignment_ref["Homo"] = Homo_desc
					assigned = 1
					break

			if not assigned:
				for temp_node_id in lineage_ref:
					temp_obj = ncbi.get_taxid_translator([temp_node_id])
					temp_name = temp_obj[temp_node_id]

					if temp_name == "Mus":
						#print("Mouse!")
						if "Mus" not in assignment_ref.keys():
							#description += "Mus\t"+hit_ref.id+"\t"
							target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
							Mus_desc = "\t".join(["Mus",hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
							assignment_ref["Mus"] = Mus_desc 
						assigned = 1
						break

			if not assigned:
				if "other" not in assignment_ref.keys():
					target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
					lineage_desc = "\t".join([lineage,hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
					assignment_ref["other"] = lineage_desc
				assigned = 1

	# check to see if it is bacteria sequence
	if len(lineage_ref) >= 2 and not assigned:
		node_id = lineage_ref[1]
		obj = ncbi.get_taxid_translator([node_id])
		name = obj[node_id]
		#print(name)
		if name == "Bacteria":
			#print("Bacteria!")
			if "Bacteria" not in assignment_ref.keys():
				target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
				lineage_desc = "\t".join([lineage,hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
				#print(lineage_desc)
				assignment_ref["Bacteria"] = lineage_desc
			assigned = 1

	# check to see if it is a phage virus sequence
	if not assigned:
		node_id = lineage_ref[0]
		obj = ncbi.get_taxid_translator([node_id])
		name = obj[node_id]
		#print(name)
		if name == "Viruses":
			#print("Virus!")
			for temp_node_id in lineage_ref:
				temp_obj = ncbi.get_taxid_translator([temp_node_id])
				#print(temp_obj)
				temp_name = temp_obj[temp_node_id]
				#print(temp_name)
				#description += temp_name+";"
				if temp_name.lower() in phage_lowercase:
					#print("phage hit!")
					if "Phage" not in assignment_ref.keys():
						#print(hit_ref.description)
						target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
						lineage_desc = "\t".join([lineage,hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
						#print(lineage_desc)
						assignment_ref["Phage"] = lineage_desc
						#print(assignment_ref["Phage"])
					assigned = 1
					break

	# check to see if it is a virus sequence
	description = ""
	if not assigned:
		print("Not assigned")
		node_id = lineage_ref[0]
		obj = ncbi.get_taxid_translator([node_id])
		name = obj[node_id]
		#print(name)
		if name == "Viruses":

			if "Viruses" not in assignment_ref.keys():
				target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
				lineage_desc = "\t".join([lineage,hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
				assignment_ref["Viruses"] = lineage_desc
				#print(assignment_ref["Viruses"])
			assigned = 1

	# check to see if it is a fungi sequence
	if len(lineage_ref) >= 4 and not assigned:
		node_id = lineage_ref[3]
		obj = ncbi.get_taxid_translator([node_id])
		name = obj[node_id]
		#print(name)
		if name == "Fungi":
			if "Fungi" not in assignment_ref.keys():
				target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
				lineage_desc = "\t".join([lineage,hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
				assignment_ref["Fungi"] = lineage_desc
			assigned = 1

	# if still not assigned, assigned to "other" category
	if not assigned:
		if "other" not in assignment_ref.keys():
			target_span = "["+str(hit_ref.hsps[0].hit_start)+"\t"+str(hit_ref.hsps[0].hit_end)+"]"
			lineage_desc = "\t".join([lineage,hit_ref.accession,str(hit_ref.seq_len),hit_ref.description,str(hit_ref.hsps[0].aln_span),str(hit_ref.hsps[0].ident_pct)+"%",target_span,str(hit_ref.hsps[0].evalue)])
			assignment_ref["other"] = lineage_desc
		assigned = 1
	return assigned,assignment_ref

###################################################################################

# Don't need to change anything below

# get file name prefix

#file_name_prefix = file_name_prefix.replace(".blastnv.out", "")
file_prefix_arr = blastout.split(".")
print(file_prefix_arr)
#file_name_prefix = file_prefix_arr[0]
file_name_prefix = ".".join(file_prefix_arr[0:3])
#file_name_suffix = "."+str(file_prefix_arr[-1])
if file_prefix_arr[3] == "mmseqsv":
	outFile = bdir+"/"+file_name_prefix+".mmseqsv.parsed"
	e_cutoff = 1e-5
	extra_rem = False
	out_suffix = ".MMseqs_VIRUSDB"
	os.system('sed -i "$ d" {0}'.format(bdir+"/"+blastout))
	vhunter_nuc = True
elif file_prefix_arr[3] == "BLASTX_VIRUSDB":
	outFile = bdir+"/"+file_name_prefix+".BLASTX_VIRUSDB.parsed"
	e_cutoff = 1e-3
	extra_rem = True
	out_suffix = ".BXVIRUSDB"
	vhunter_nuc = False

# create ouput file

try:
    out = open(outFile, 'w')
except IOError:
	print("can not open file "+outFile)

# create a tmp_taxonomy directory in the x directory if tmp_taxonomy does not exist
#os.system('mkdir -p ~/mmseqs_TESTING/databases/tmp_taxonomy')


#with open("/Users/qthomas/mmseqs_TESTING/database/phage_list.txt", "r") as phageFile:
#	phage_list = phageFile.read().splitlines()

phage_list = ["Uroviricota", "Loebvirae","Trapavirae","Sangervirae","Caudovirales","Caudoviricetes","unclassified bacterial viruses"]
phage_lowercase=[x.lower() for x in phage_list]

#print(phage_lowercase)

unassigned = [] # query should be kept for further analysis
virusseq = []
phageseq = []
total_records = 0

print("parsing blast output files...\n\n")

input_file = bdir+"/"+blastout
custom_fields=["qseqid","sacc","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sseqid","qlen","slen"]
report = SearchIO.parse(input_file, "blast-tab", fields=custom_fields)
#print(report)

# Go through BLAST reports one by one
for result in report:
	#print(result.__dir__())
	total_records+=1
	haveHit = 0
	keep_for_next_step = 1
	assignment = {}
	assignment_NotBestE = {}
	# only take the best hits
	best_e = 100
	hit_count = 0
	determined = 0
	#print(hit)

	for hit in result:
		hit_desc = hit.id
		#print(hit_desc)
		hit.description = hit_desc.split(" ", 1)[1]
		#print(hit.__dir__())
		# from hit name get hit gi number
		hit_name = hit.accession
		#print(hit.id)
		#print(hit.description)
		haveHit = 1
		hit_count+=1
		if hit_count == 1:
			best_e = hit[0].evalue
			#print(best_e)

		# check whether the hit should be kept for further analysis
		if best_e <= e_cutoff: # similar to known, need Phylotyped
			#print(hit[0].evalue)
			if hit[0].evalue == best_e: # only get best hits
				#get taxonomy lineage
				#print(hit_name)
				if vhunter_nuc:
					sth = cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_nucl where accession_version = '"+hit_name+"'")
				else:
					sth = cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_prot where accession_version = '"+hit_name+"'")
				ref = sth.fetchone()
				#print(ref)
				if ref: # some gi don't have record in gi_taxid_nucl
					taxID = ref[2]
					taxon_name = ncbi.get_taxid_translator([taxID])
					# print(taxon_name)
					if not taxon_name:
						description = "undefined taxon "+result.description+"\t"+hit.id+"\t"+hit[0].evalue
						assignment["other"] = description
					else:
						lineage = ncbi.get_lineage(taxID)
						lineage.pop(0)
						if lineage:
							print(lineage)
							determined = 1
							success,assignment = PhyloType(lineage, result, hit)

				else: # for situations that gi does not have corresponding taxid
					determined = 1
					target_span = "["+str(hit.hsps[0].hit_start)+"\t"+str(hit.hsps[0].hit_end)+"]"
					result_desc = "\t".join([hit.description,hit.accession,str(hit.seq_len),hit.description,str(hit.hsps[0].aln_span),str(hit.hsps[0].ident_pct)+"%",target_span,str(hit.hsps[0].evalue)])
					assignment["other"] = result_desc

			elif hit[0].evalue <= e_cutoff: # significant but is not the same e value as best hit
				# get taxonomy lineage
				if vhunter_nuc:
					sth = cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_nucl where accession_version = '"+hit_name+"'")
				else:
                    			sth = cursor_dbhsqlite.execute("SELECT * FROM acc_taxid_prot where accession_version = '"+hit_name+"'")
				ref = sth.fetchone()
				if ref: # some gi don't have record in gi_taxid_nucl
					taxID = ref[2]
					taxon_name = ncbi.get_taxid_translator([taxID])

					####################################
					#I have not actually seen an example of this in any of the outputs so I'm not sure how to test it?
					if not taxon_name:
						target_span = "["+str(hit.hsps[0].hit_start)+"\t"+str(hit.hsps[0].hit_end)+"]"
						undef_desc = "\t".join(["undefined taxon ",hit.accession,str(hit.seq_len),hit.description,str(hit.hsps[0].aln_span),str(hit.hsps[0].ident_pct)+"%",target_span,str(hit.hsps[0].evalue)])
						assignment_NotBestE["other"] = description
					####################################
					else:
						lineage = ncbi.get_lineage(taxID)
						lineage.pop(0)
						if lineage:
							success,assignment_NotBestE = PhyloType(lineage, result, hit)

						#################################################################
						# If the sequence also hit any other species with significant e value skip all the rest hits.
						#if (((defined $assignment_NotBestE{"Bacteria"}) || (defined $assignment_NotBestE{"Fungi"}) || (defined $assignment_NotBestE{"Homo"}) || (defined $assignment_NotBestE{"Mus"}) || (defined $assignment_NotBestE{"Phage"}) || (defined $assignment_NotBestE{"other"})) ) {
						if any(x in ["Bacteria", "Fungi", "Homo", "Mus", "Phage", "other"] for x in assignment_NotBestE.keys()):
							break

		# finish phylotype for given hit
		else: # e value is not significant
			if determined: # skip the rest hits that are not significant
				break
			else:
				target_span = "["+str(hit.hsps[0].hit_start)+"\t"+str(hit.hsps[0].hit_end)+"]"
				nosig_desc = "\t".join(["hit not significant",hit.accession,str(hit.seq_len),hit.description,str(hit.hsps[0].aln_span),str(hit.hsps[0].ident_pct)+"%",target_span,str(hit.hsps[0].evalue)])
				assignment["unassigned"] = nosig_desc
				break
	# finish all hits
	if not haveHit:
		assignment["unassigned"] = "no hit"

	
	# remove duplicate assignment
	# If a query is assigned both Homo and Primates, it will be reported as Homo only
	# If a query is assigned a real taxon name and "other" for reason like"other sequences;
	# artificial sequences", or no taxon id in taxon database it will be reported only as 
	# the real taxon name
	#print(assignment.keys())
	num_assignment = assignment.keys()
	if len(num_assignment) > 1: # have multiple assignment, can only be virus or phage
		###############################################
		# determine phage hits
		# If a sequence hits virus and  phage, the sequence is assigned to "Phage" category. 
		if extra_rem:
			if "Viruses" in assignment.keys():
				if "Phage" in assignment.keys():
					del assignment["Viruses"]
				if "other" in assignment.keys():
					del assignment["other"]
				if "unassigned" in assignment.keys():
					del assignment["unassigned"]
			if "Phage" in assignment.keys():
				if "other" in assignment.keys():
					del assignment["other"]
				if "unassigned" in assignment.keys():
					del assignment["unassigned"]
		else:
			if "Viruses" in assignment.keys():
				if "Phage" in assignment.keys():
					del assignment["Viruses"]
	
	elif len(num_assignment) == 1: # have exactly one assignment
		if "Viruses" in assignment.keys(): # it's virus assignment
			if "Phage" in assignment_NotBestE.keys():# but has phage as significant (not best) hit
				#print("has phage hits!!!!!!!!!!!!!!\n")
				assignment["Phage"] = assignment_NotBestE["Phage"] # considered to be phage seq
				del assignment["Viruses"]
	assign_to_virus=0
	for assign in assignment.keys():
		if assign == "Viruses":
			virusseq.append(result.id)
			assign_to_virus=1
		if assign == "Phage":
			phageseq.append(result.id)
			assign_to_virus=1
		out.write(result.id+"\t"+str(result.seq_len)+"\t"+assign+"\t"+assignment[assign]+"\n")

	if assign_to_virus==0:
		unassigned.append(result.id)



if total_records == 0:
        percent_unassigned = 0.0
else:
        percent_unassigned = (len(unassigned)*100)/total_records

out.write("# Summary: "+ str(len(unassigned))+" out of "+str(total_records)+" ("+str(format(percent_unassigned, ".2f"))+"%) are unassigned.\n")
out.close()

file = bdir+"/"+file_name_prefix+".fasta"
seq_dict = read_FASTA_data(file)

outFile = bdir+"/"+file_name_prefix+out_suffix+"_hit.fa"
try:
	out2 = open(outFile, "w")
except:
	print("Can't create output file : "+outFile)
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
for seq_name in unassigned:
	out2.write(">"+seq_name+"\n")
	seq=seq_dict[seq_name].seq
	out2.write(str(seq)+"\n")
out2.close()

