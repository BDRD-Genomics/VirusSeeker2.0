#!/usr/bin/env perl
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

use strict;
use warnings;
use Config::Simple;

my $version = 0.04;
 

# Get Paths from VS.config file
my $config_path = `dirname $0`;
chomp $config_path;
$config_path = $config_path."/VS.config";
my $cfg = new Config::Simple($config_path);

# SLURM settings
my $partition = $cfg->param('slurm_partition');

# path and name of databases
my $db_MMSEQS_NT = $cfg->param('mmseqsnt');
my $mmseqs_tmpdir = $cfg->param('mmseqs_tmpdir');
my $db_NR = $cfg->param('nr');
my $db_MMSEQS_NT_VIRUS = $cfg->param('virus_nt');
my $db_NR_VIRUS = $cfg->param('virus_nr');

# bin path for any binaries not in conda env (currently only RepeatMasker)
my $bin_path = $cfg->param('bin_path');

# my reference database path
my $ref_path = $cfg->param('ref');
my $viral_reference_genome = $ref_path."/ref_viruses_rep_genomes";
 
# my conda prefix
my $conda_prefix = $cfg->param('conda_prefix');

#color code
my $red = "\e[31m";
my $green = "\e[32m";
my $yellow = "\e[33m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $gray = "\e[37m";
my $normal = "\e[0m"; 

# usage information
(my $usage = <<OUT) =~ s/\t+//g;
This script will run the VirusSeeker pipeline on SLURM.
Pipeline version: $version
$yellow		Usage: perl $0 <sample_folder> <ref genome> <#CPU> <Memory(GB)> <use checkpointing> <step_number> <VStrack> $normal

<sample_folder> = full path of the folder holding files for the sample
<ref genome> = 1. Human genome
 
<#CPU> = number of CPU to use
<Memory> = GB of memory to use
<use checkpointing> = 1. yes, 0. no
<step_number>** = [1..39] run this pipeline step by step. $red **If <step_number> is set to -1, remove adapter and host removal skipped; running the whole pipeline if <step_number> is 0.
<VStrack>** = run Discovery or Virome pipeline. $red **If <VStrack> is set to V, assembly is skipped; run the whole pipeline if <VStrack> is D.

$yellow	[1]   Remove Adapter 
		[2]   Stitch SE1 and SE2
		[3]   Sequence quality control 

$purple If host removal is ON, [4, 5] are:
		[4]   Align to reference genome ( BWA-MEM )
		[5]   Get unmappped reads

$cyan If host removal is turned OFF [4,5] are: 
		[4]   fastq to fasta and renaming of files to submit to CD-HIT
		[5]   no job 5 if host removal is turned off
		
$yellow	[6]   run CD-HIT to cluster
		[7]   Extract the top 3 similar reads 

		[7]   Extract original reads to FASTQ format
		[9]   Split for assembly 
		[10]  Assembly with Newbler 
		[11]  Extract Contig, Singleton and Outlier from assembly
		[12]  Split above sequences for a second round assembly
		[13]  Second round assembly, contigs use Phrap, Singleton and Outlier sequences use Newbler 
		[14]  Second Extraction of Singleton and Outlier from the previous step assembly

		[15]  Get all contigs, Singletons, Outliers
		[16]  Run cd-hit on all contigs, Singletons, Outliers

		[17]  Prepare for RepeatMasker
		[18]  Submit RepeatMasker job array
		[19]  Sequence quality control after RepeatMasker

		[20] Split files for MMseqs Virus-Only Database
		[21] Submit MMseqs Virus-Only Database job array
		[22] Parse Viral MMseqs result

		[23] Split files for Diamond BlastX Virus-Only Database
		[24] Submit Diamond BlastX Virus-Only Database job array
		[25] Parse Viral Diamond BlastX result

		[26] Split for MMseqs against NT 
		[27] Run MMseqs search
		[28] Parse MMseqs output file 

		[29]  Split for MMseqs NT
		[30]  Submit MMseqs job array
		[31]  Parse MMseqs result

		[32]  Pool and split files for Diamond BlastX
		[33]  Submit Diamond BlastX job array
		[34]  Parse Diamond BlastX result

		[35] Generate assignment report for the sample
		[36] Generate assignment summary for the sample

		[37] Generate phage report
		[38] Generate phage summary
		[39] Run VS Read Counter

$red To generate tabular results, must run Generate_assignment_summary.pl - use sample directory as input
		
$normal
OUT

die $usage unless scalar @ARGV == 9;
my ($sample_dir, $ref_genome_choice, $num_CPU, $slurm_mem, $assembly_type, $assembly_mode, $use_checkpoint, $step_number, $VStrack) = @ARGV;
die $usage unless (($step_number >= -1)&&($step_number <= 32));
die $usage unless (($VStrack eq "D")||($VStrack eq "V"));
if ($num_CPU == 0) {
    # Hardware resource defaults
    $num_CPU = $cfg->param('slurm_cpus_per_task');
    $slurm_mem = $cfg->param('slurm_mem');
}
my $num_input_files = 0;
if ($assembly_type eq "s") {
	$num_input_files = 2;
}
if ($assembly_type eq "l") {
	$num_input_files = 1;
}
if ($assembly_type eq "h") {
	$num_input_files = 3;
	$num_CPU = "64";
	$slurm_mem = "128g";
}
print "assembly type=$assembly_type\n";
print "num_input_files=$num_input_files\n";
my reference_genome=$ref_genome_choice

#my $arg_log = "arguments".$sample_name."_out.txt";
#open(STCH, ">$sample_dir/$arg_log") or die $!;
#print STCH "sample_dir = $sample_dir \n ref_genome = $ref_genome_choice \n num_CPU = $num_CPU \n assembly = $assembly_type";  
#close STCH;

#####################################################################################
# get name of the sample and path to the data
if ($sample_dir =~/(.+)\/$/) {
	$sample_dir = $1;
}
my $sample_name = (split(/\//,$sample_dir))[-1];
my $account = "VS:".(split(/\//,$sample_dir))[-2];
my $sample_path = $account."/".$sample_name;
print $sample_name,"\n";
my $last_jobid1 ="";
my $last_jobid2 ="";

my $arg_log = "arguments_".$sample_name."_out.txt";
open(STCH, ">$sample_dir/$arg_log") or die $!;
print STCH "sample_dir = $sample_dir \n ref_genome = $ref_genome_choice \n num_CPU = $num_CPU \n assembly = $assembly_type\n mode = $VStrack\n num_input_files = $num_input_files\n assembly_mode = $assembly_mode\n\n";
close STCH;


####################################################################################
# Here are the parameters that can be adjusted according to computer cluster 
# configuration.
# Use this as max number of reads per file as input for assembler. 
####################################################################################
# To run jobs faster, split large fasta files to small ones. Split to specific number 
# of files instead of specific sequences in each small file, because the number of job 
# array cannot be determined if spliting to specific number of sequences in each file. 
# Job number is required by qsub ${SLURM_ARRAY_TASK_ID} the minimum size of each file is 4kb

# Repeatmasker
# It takes ~ 0.4 sec to process each seq. 1200 seq takes ~ 8 min.
# minimum number of sequences per file when spliting into smaller files
my $num_seq_per_file_Repeatmasker = 36000;

# the number of job files for RepeatMasker. The most optimal number would be an 
# estimated number that makes the number of sequences in each file close to 2400. 
#my $file_number_of_RepeatMasker = 1000; 
my $file_number_of_RepeatMasker = 500; 

# MMseqs2 search
# 20,000 reads use ~15 min in MMseqs2 against reference genome. If too many reads in 
# one file the memory usage is too big. Job may be killed before it finishes. 
#my $num_seq_per_file_MMseqs_host = 40000;
# If 1 sample is sequenced in 1 MiSeq run, 400 file would be max number of files
# based on percentage of reads left for analysis at each step from MiSeq run 1 data. 
# This number should be bigger to make sure each file does not have too many reads.


my $contig_length_cutoff = 500; # length cutoff for a contig to be in output 

# MMseqs virus only database. 10,000 sequence takes ~ 10 min
# Specify the number of sequences to be in each file 
my $num_seq_per_file_MMseqs_VIRUSDB = 5000000;
# Specify the total number of files to be splited into
my $file_number_of_MMseqs_VIRUSDB = 50; 

# Diamond BLASTX virus only database.
# Specify the number of sequences to be in each file 
my $num_seq_per_file_BLASTX_VIRUSDB = 5000;
# Specify the total number of files to be splited into
my $file_number_of_BLASTX_VIRUSDB = 50; 

# MegaBLAST NT
# 2,000 reads use ~40 min in MegaBLAST against reference genome. If too many reads in 
# one file the memory usage is too big. Job may be killed before it finishes. 
my $num_seq_per_file_MMseqs_NT = 50000000;
# If 1 sample is sequenced in 1 MiSeq run, 400 file would be max number of files
# based on percentage of reads left for analysis at each step from MiSeq run 1 data. 
# This number should be bigger to make sure each file does not have too many reads.
my $file_number_of_MMseqs_NT = 25;


# BLASTN: 1000 reads takes ~ 6 hours. 
# minimum number of sequences per file
my $num_seq_per_file_MMseqs = 500000000;
# Specify the number of small fasta files to split to from a large file for BLASTN
# The most optimal number would be an estimated number that makes the number of 
# sequences in each file close to 200.
#my $file_number_of_BLASTN = 1200; 
my $file_number_of_MMseqs = 25; 

# BLASTX: 1000 reads takes ~ 15 hours.
# Specify the number of sequences to be in each file for BLASTX NR
my $num_seq_per_file_BLASTX = 50;
# Specify the number of small fasta files to split to from a large file for BLASTX
#my $file_number_of_BLASTX = 800; 
my $file_number_of_BLASTX = 50; 

####################################################################################
# software paths
my $mmseqs_easylinclust = "mmseqs easy-linclust";
my $mmseqs_easysearch = "mmseqs easy-search";
my $mmseqs_touchdb = "mmseqs touchdb";
my $mmseqs_cutoff = 0.95;
my $mmseqs_parameter = "--min-seq-id $mmseqs_cutoff -c $mmseqs_cutoff --threads $num_CPU";
my $repeat_masker = $bin_path."RepeatMasker/RepeatMasker";
my $blastn = "export BLAST_USAGE_REPORT=false; blastn";
my $diamond_args = $cfg->param('diamond_args');
my $diamond = "diamond blastx $diamond_args";
my $bbduk = "bbduk.sh";
my $samtools = "samtools";
my $bamToFastq = "bamToFastq";
my $cutadapt = "cutadapt";
my $pear = "pear";
my $seqtk = "seqtk";
my $pigz = "pigz";
my $reformat = "reformat.sh";
my $SPAdes_mem = $slurm_mem =~ s/g//gr;
if ($SPAdes_mem < 250) { $SPAdes_mem = 250; }
my $metaSPAdes = "spades.py --only-assembler --tmp-dir /dev/shm --meta";
my $SPAdes = "spades.py --only-assembler --tmp-dir /dev/shm --isolate";
my $seqkit = "seqkit grep -v -f";
my $dragonflye = "dragonflye --force";
my $unicycler = "unicycler";
my $minimap = "minimap2 -ax map-ont";

####################################################################################
# folders to be created to store intermediate results
my $MMseqs_dir_host1 = $sample_dir."/".$sample_name."_SE1_MMseqs_HOST";
my $MMseqs_dir_host2 = $sample_dir."/".$sample_name."_SE2_MMseqs_HOST";
my $assembly_dir_metaSPAdes = $sample_dir."/".$sample_name."_metaSPAdes";
my $assembly_dir_dragonflye = $sample_dir."/".$sample_name."_dragonflye";
my $assembly_dir_unicycler = $sample_dir."/".$sample_name."_unicycler";
my $seq_dir = $sample_dir."/".$sample_name."_all_seqs";
my $repeatmasker_dir = $sample_dir."/".$sample_name."_RepeatMasker";
my $MMseqs_VIRUSDB_DIR_SUFFIX = "_MMseqs_VIRUSDB"; # blastn against virus-only (NT extracted) database
my $BLASTX_VIRUSDB_DIR_SUFFIX = "_BLASTX_VIRUSDB"; # blastx against virus-only (NR extracted) database
my $MMseqs_dir_NT = $sample_dir."/".$sample_name."_MMseqs_NT";
my $MMseqs_NT_dir = $sample_dir."/".$sample_name."_MMseqs_NT";
my $BLASTX_NR_dir =$sample_dir."/".$sample_name."_BLASTX_NR";

####################################################################################
# Everything else below should be automated.
my $HOME = $ENV{HOME};

# store job files here
my $job_files_dir = $sample_dir."/job_script";
if (! -d $job_files_dir) {
	`mkdir -p $job_files_dir`;
}

#create a folder to store SGE output and error files
my $sge_files_dir = $sample_dir."/SGE_DIR";
#if (! -d $sge_files_dir) {
	#`mkdir -p $sge_files_dir`;
#}

#create a folder to store SLURM output and error files
my $SLURM_files_dir = $sample_dir."/SLURM_DIR";
if (! -d $SLURM_files_dir) {
	`mkdir -p $SLURM_files_dir`;
}

#create a folder to store status files
my $status_log = $sample_dir."/status_log";
if (! -d $status_log) {
	`mkdir -p $status_log`;
}

my $run_script_path = `dirname $0`;
my $bash_bin_path = `dirname $0`;
chomp $bash_bin_path;
chomp $run_script_path;
$run_script_path = "perl ".$run_script_path."/";
my $qsub_com = "";

my $current_job_file1 = "no";#cannot be empty
my $current_job_file2 = "no";#cannot be empty
my $hold_job_file1 = "";
my $hold_job_file2 = "";

####################################################################################
# run the whole pipeline

#short + meta = metaspades
#short + iso = spades
#long + meta = dragonflye
#long + iso = dragonflye
#hybrid + meta = dragonflye
#hybrid + iso = dragonflye + unicycler

if ($step_number == 0 or $step_number == -1) {
	# 1
	#&remove_adapter();
        # 2
        &quality_control();
        # 3
        if ($step_number == 0) {
        	if ($assembly_type eq "s" or $assembly_type eq "h"){
                &map_to_host();
            }
            if ($assembly_type eq "l" or $assembly_type eq "h"){
                &map_to_host_minimap();
            }
        }elsif ($step_number == -1) {
                &skip_host_removal();
        }
        
        
        if ($VStrack eq "V") {
                # 4
                &skip_assembly();
        }else {
		
		if ($assembly_type eq "s" and $assembly_mode eq "m"){
				# 4
                		&metaSPAdes_assembly();
                		# 5
                		&map_to_assembly();
                }
                if ($assembly_type eq "s" and $assembly_mode eq "i"){
                                # 4
                                &SPAdes_assembly();
                                # 5
                                &map_to_assembly();
                }
                if ($assembly_type eq "l"){
						# 4
                		&dragonflye_assembly();
                		# 5
                		&map_to_assembly_minimap();
                }
                if ($assembly_type eq "h" and $assembly_mode eq "i"){
						# 4
                		&unicycler_assembly();
				&dragonflye_assembly();
                		# 5
                		&map_to_assembly();
                		&map_to_assembly_minimap();
                }
                if ($assembly_type eq "h" and $assembly_mode eq "m"){
                                                # 4
                                &dragonflye_assembly();
                                # 5
                                &map_to_assembly();
                                &map_to_assembly_minimap();
                }
        }
        # 6
        if ($assembly_type eq "s" or $assembly_type eq "h"){
        	&stitching();
        }
        # 7
        &get_all_seqs();
        # 8
        &run_mmseqs();
        # 9
        &pre_RepeatMasker();
        # 10
        &RepeatMasker();
        # 11
        &post_RepeatMasker();
        # 12
        &split_for_MMseqs_VIRUSDB();
        # 13
        &submit_job_array_MMseqs_VIRUSDB();
        # 14
        &parse_MMseqs_VIRUSDB();
        # 15
        &split_for_BLASTX_VIRUSDB();
        # 16
        &submit_job_array_BLASTX_VIRUSDB();
        # 17
        &parse_BLASTX_VIRUSDB();
        # 18
        &split_for_MMseqs_NT();
        # 19
        &submit_job_array_MMseqs_NT();
        # 20
        &parse_MMseqs_NT();
        # 21
        #&pool_split_for_MMseqs();
        # 22
        #&submit_job_array_MMseqs();
        # 23
        #&parse_MMseqs();
        # 21
        &pool_split_for_BLASTX();
        # 22
        &submit_job_array_BLASTX();
        # 23
        &parse_BLASTX();
        # 24
	if ($assembly_type eq "s" or $assembly_type eq "h"){
	        &map_to_viral_ref();
        }
	if ($assembly_type eq "l" or $assembly_type eq "h"){
                &map_to_viral_ref_b();
        }

	# 25
        &generate_assignment_report();
        # 26
        &generate_assignment_summary();
        # 27 
        &generate_phage_report();
        # 28
        &generate_phage_summary();
        # 29
        &generate_supplemental_read_counts();

}elsif ($step_number == 1) {
	&remove_adapter();
}elsif ($step_number == 2) {
        &quality_control();
}elsif ($step_number == 3) {
        &map_to_host();
}elsif ($step_number == 4) {
        &metaSPAdes_assembly();
}elsif ($step_number == 5) {
        &map_to_assembly();
}elsif ($step_number == 6) {
        &stitching();
}elsif ($step_number == 7) {
        &get_all_seqs();
}elsif ($step_number == 8) {
        &run_mmseqs();
}elsif ($step_number == 9) {
        &pre_RepeatMasker();
}elsif ($step_number == 10) {
        &RepeatMasker();
}elsif ($step_number == 11) {
        &post_RepeatMasker();
}elsif ($step_number == 12) {
        &split_for_MMseqs_VIRUSDB();
}elsif ($step_number == 13) {
        &submit_job_array_MMseqs_VIRUSDB();
}elsif ($step_number == 14) {
        &parse_MMseqs_VIRUSDB();
}elsif ($step_number == 15) {
        &split_for_MMseqs_VIRUSDB();
}elsif ($step_number == 16) {
        &submit_job_array_BLASTX_VIRUSDB();
}elsif ($step_number == 17) {
        &parse_BLASTX_VIRUSDB();
}elsif ($step_number == 18) {
        &split_for_MMseqs_NT();
}elsif ($step_number == 19) {
        &submit_job_array_MMseqs_NT();
}elsif ($step_number == 20) {
        &parse_MMseqs_NT();
}
#elsif ($step_number == 21) {
#        &pool_split_for_MMseqs();
#}
#elsif ($step_number == 22) {
#        &submit_job_array_MMseqs();
#}elsif ($step_number == 23){
#        &parse_MMseqs();
elsif ($step_number == 21) {
        &pool_split_for_BLASTX();
}elsif ($step_number == 22) {
        &submit_job_array_BLASTX();
}elsif ($step_number == 23) {
        &parse_BLASTX();
}elsif ($step_number == 24) {
        &map_to_viral_ref();
}elsif ($step_number == 25){
        &generate_assignment_report();
}elsif ($step_number == 26){
        &generate_assignment_summary();
}elsif ($step_number == 27) {
        &generate_phage_report();
}elsif ($step_number == 28) {
        &generate_phage_summary();
}elsif ($step_number == 29) {
        &generate_supplemental_read_counts();
}
else{
        #die $VStrack;
        die $usage;
}

exit;

##########STAGE 1: Remove adapter sequences from reads##########
sub remove_adapter {
	$current_job_file1 = "j1_".$sample_name."_RemoveAdapter1.sh";


	my $Adapter_file = $sample_dir."/".$sample_name."_adapter.txt";
	my @adapter_seqs = ();
	open(IN, "<$Adapter_file") or die "$Adapter_file does not exist!!!";
	while (<IN>) {
		if ($_ =~ /^\s$/) { # skip blank line
			next;
		}
		chomp;
		push @adapter_seqs, $_;
	}
	close IN;

	my $command_line_adapter_seq = "";
	foreach my $adapter (@adapter_seqs) {
		$command_line_adapter_seq .= " -b \"".$adapter."\" ";
	}
	 
	######################################################################
	open(RemoveAdapter, ">$job_files_dir/$current_job_file1") or die $!;
	print RemoveAdapter "#!/bin/bash\n";
	print RemoveAdapter "#SBATCH --partition=$partition\n";
	print RemoveAdapter "#SBATCH --account=$account\n";
	print RemoveAdapter "#SBATCH --job-name=${sample_path}:stage1\n";
	print RemoveAdapter "#SBATCH --array=1-2\n";
	print RemoveAdapter "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	print RemoveAdapter "IN=$sample_dir"."/".$sample_name."_SE\${SLURM_ARRAY_TASK_ID}.fastq.gz\n";
	print RemoveAdapter "OUTFILE=$sample_dir"."/".$sample_name."_SE\${SLURM_ARRAY_TASK_ID}.RemoveAdapter.fastq.gz\n";
	print RemoveAdapter "REPORT=$sample_dir"."/".$sample_name."_SE\${SLURM_ARRAY_TASK_ID}.RemoveAdapter.report\n\n";
        print RemoveAdapter "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print RemoveAdapter "conda activate vs\n";  

	# Check if RemoveAdapter finished successfully. If true, resubmitting will skip this step and keep results
	print RemoveAdapter "if [ ! -e $status_log/j1_RemoveAdapter\${SLURM_ARRAY_TASK_ID}_finished ]\n";
	print RemoveAdapter "then\n";
	print RemoveAdapter " 	$cutadapt ".$command_line_adapter_seq." -n 5 -O 5  -o \${OUTFILE}  \${IN} > \${REPORT}   \n";
	print RemoveAdapter "	OUT=\$?\n"; 
	print RemoveAdapter '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print RemoveAdapter "	then\n";
	print RemoveAdapter "		echo \"Fatal Error trying to run cutadapt.\"  \n";
	print RemoveAdapter "		exit 1\n";
 
	# if prinseq finished succussfully, make a file named QC_finished
	print RemoveAdapter "	else\n";

	# RemoveAdapter finished succussfully, make a file named RemoveAdapter1_finished
	print RemoveAdapter "		touch $status_log/j1_RemoveAdapter\${SLURM_ARRAY_TASK_ID}_finished\n";
	print RemoveAdapter "	fi\n\n";
	print RemoveAdapter "fi\n\n";
	close RemoveAdapter;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 2: Quality control with bbduk##########
sub quality_control{
        $current_job_file1 = "j2_".$sample_name."_QC.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage2\n";
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";

        #if (!$step_number or $step_number == -1) {
        #        print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        #}
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  

        #print STCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.RefGenome.unmapped.fastq.gz\n";
        print STCH "IN1=$sample_dir"."/".$sample_name."_SE1.fastq.gz\n";
        print STCH "IN2=$sample_dir"."/".$sample_name."_SE2.fastq.gz\n";
        print STCH "IN3=$sample_dir"."/".$sample_name.".fastq.gz\n";
        print STCH "OUTFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.fastq.gz\n";
        print STCH "OUTFILE2=$sample_dir"."/".$sample_name.".RemoveAdapter.LR.fastp.fastq.gz\n";
        print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".bbduk_time.txt\n";

        print STCH "mkdir -p $seq_dir \n\n";
        print STCH "if [ ! -e $status_log/j2_QC_finished ]\n";
        print STCH "then\n";
        print STCH "date > \${TIMEFILE}\n";
        
        print STCH "	if [ $num_input_files -eq 2 ] || [ $num_input_files -eq 3 ]\n";
        print STCH "	then\n";
        my $com1 = $bbduk." in1=\${IN1} in2=\${IN2} out=\${OUTFILE} -Xmx${slurm_mem} tossbrokenreads ref=phix,adapters tbo tpe qtrim=rl trimq=30 maq=20 minlength=50 maxns=6 minbasequality=10";
        print STCH "    	", $com1, "\n";
        print STCH "	fi\n";
        print STCH "	if [ $num_input_files -eq 3 ] || [ $num_input_files -eq 1 ]\n";
        print STCH "	then\n";
		print STCH "		fastp --dedup --overrepresentation_analysis --low_complexity_filter --trim_poly_x --report_title '${sample_name}_fastp_report' -q 10 -e 10 --cut_front --cut_tail -w 16 -j ${sample_dir}/${sample_name}_fastp.json -h ${sample_dir}/${sample_name}_fastp.html --in1 \${IN3} -o \${OUTFILE2}\n";
        print STCH "    fi\n";

        # Check if bbduk finished successfully. If true, resubmitting will skip this step and keep results
        print STCH "    OUT=\$?\n";
        print STCH '    if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
        print STCH "    then\n";
        print STCH "            echo \"Fatal Error trying to run bbduk.\"  \n";
        print STCH "            exit 1\n";

        # If bbduk finished succussfully, make a file named QC_finished
        print STCH "    else\n";
        print STCH "            date >> \${TIMEFILE}\n";
        print STCH '            echo "', $com1, '"  >> ${TIMEFILE} ', "\n";
        print STCH "            touch $status_log/j2_QC_finished\n";
        print STCH "    fi\n";
        print STCH "fi\n\n";
        close STCH;

        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 3a: map reads to reference genome##########
sub map_to_host{
        # this is the job script
        my $job_script = $job_files_dir."/j3_".$sample_name."_map_to_host_script.sh";
        open(JSTCH, ">$job_script") or die $!;
        chmod 0755, $job_script;
        print JSTCH "#!/bin/bash\n";
        print JSTCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.fastq.gz\n";
        print JSTCH "OUTFILE1=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.mapped.fastq.gz\n";
        print JSTCH "OUTFILE2=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.fastq.gz\n";
        print JSTCH "OUTFILE3=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n";
        print JSTCH "OUTFILE4=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n";
        print JSTCH "OUTFILE5=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n";
        print JSTCH "TIMEFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.map_time.txt\n";
        print JSTCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print JSTCH "conda activate vs\n";

        # Check if mapping finished successfully. If true, resubmitting will skip this step and keep results
        print JSTCH "if [ ! -e $status_log/j3_map_to_host_finished ]\n";
        print JSTCH "then\n";
        print JSTCH "   date > \${TIMEFILE}\n";
        print JSTCH "   bbmap.sh in=\${IN} ref=$reference_genome tossbrokenreads covstats=$sample_dir"."/".$sample_name.".bbmap_host_covstats.txt statsfile=$sample_dir"."/".$sample_name.".bbmap_host_statsfile.txt -Xmx${slurm_mem} fast=t ow 32bit=t pigz=t outm=\${OUTFILE1} outu=\${OUTFILE2} \n";
        print JSTCH "   OUT=\$?\n";
        print JSTCH '   if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
        print JSTCH "   then\n";
        print JSTCH "           echo \"Fatal Error trying to run bbmap.\"  \n";
        print JSTCH "           exit 1\n";
        print JSTCH "   else\n";
        print JSTCH "           $reformat in=\${OUTFILE2} int=t out1=\${OUTFILE3} out2=\${OUTFILE4} outs=\${OUTFILE5}\n";
        print JSTCH "           touch $status_log/j3_map_to_host_finished\n";
        print JSTCH "   fi\n";

        #if mapping finished succussfully, make a file named finished
        print JSTCH "   date >> \${TIMEFILE}\n";
        print JSTCH "fi\n\n";
        close JSTCH;

        # here to submit
        $current_job_file1 = "j3_".$sample_name."_map_to_host.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage3a\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n"; 
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
        if (!$step_number or $step_number == -1) {
                print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        }
        if ($use_checkpoint) {
                print STCH "checkpoint bash  $job_script \n";
        }
        else {
                print STCH " bash  $job_script \n";
        }

        close STCH;

    $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 3b: Skip host removal, send output to QC##########
sub skip_host_removal{
        $current_job_file1 = "j3_".$sample_name."_skip_host_removal.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage3\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";

        if (!$step_number or $step_number == -1) {
                print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        }
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  
        
        print STCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.fastq.gz\n";
        print STCH "OUTFILE1=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n";
        print STCH "OUTFILE2=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n";
        print STCH "OUTFILE3=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n";
	
	#Check if host removal finished successfully. If true, resubmitting will skip this step and keep results.
        print STCH "OUT=1\n";
        print STCH "if [ ! -e $status_log/j3_skip_host_removal_finished ]\n";
        print STCH "then\n";
        print STCH "    $reformat in=\${IN} int=t out1=\${OUTFILE1} out2=\${OUTFILE2} outs=\${OUTFILE3}\n";
        print STCH "    OUT=\$?\n";
        print STCH '    if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
        print STCH "    then\n";
        print STCH "            echo \"Fatal Error trying to skip host removal.\"  \n";
        print STCH "            exit 1\n";

        # if finished, make a file named _finished
        print STCH "    else\n";
        print STCH "            touch $status_log/j3_skip_host_removal_finished\n";
        print STCH "    fi\n";
        print STCH "fi\n\n";
        close STCH;

        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
 }

##########STAGE 3c: map ONT reads to reference genome##########
sub map_to_host_minimap{
        # this is the job script
        my $job_script = $job_files_dir."/j3c_".$sample_name."_map_to_host_script.sh";
        open(JSTCH, ">$job_script") or die $!;
        chmod 0755, $job_script;
        print JSTCH "#!/bin/bash\n";
        print JSTCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.LR.fastp.fastq.gz\n";
        print JSTCH "OUTFILE1=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.minimap.sam\n";
        print JSTCH "OUTFILE2=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.mapped.LR.fastq\n";
        print JSTCH "OUTFILE3=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.LR.fastq\n";
        print JSTCH "TIMEFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.LR.map_time.txt\n";
        print JSTCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print JSTCH "conda activate vs\n";

        # Check if mapping finished successfully. If true, resubmitting will skip this step and keep results
        print JSTCH "if [ ! -e $status_log/j3c_map_to_host_finished ]\n";
        print JSTCH "then\n";
        print JSTCH "   date > \${TIMEFILE}\n";
        print JSTCH "   minimap2 -ax map-ont $reference_genome \$IN > \$OUTFILE1\n";
        print JSTCH "   OUT=\$?\n";
        print JSTCH '   if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
        print JSTCH "   then\n";
        print JSTCH "           echo \"Fatal Error trying to run minimap.\"  \n";
        print JSTCH "           exit 1\n";
        print JSTCH "   else\n";
        print JSTCH "			samtools fastq -F 4 -c 6 \$OUTFILE1 > \$OUTFILE2\n";
		print JSTCH "			samtools fastq -f 4 \$OUTFILE1 > \$OUTFILE3\n";
        print JSTCH "           touch $status_log/j3c_map_to_host_finished\n";
        print JSTCH "   fi\n";

        #if mapping finished succussfully, make a file named finished
        print JSTCH "   date >> \${TIMEFILE}\n";
        print JSTCH "fi\n\n";
        close JSTCH;

        # here to submit
        $current_job_file1 = "j3c_".$sample_name."_map_to_host.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage3c\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n"; 
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
        if (!$step_number or $step_number == -1) {
                print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        }
        if ($use_checkpoint) {
                print STCH "checkpoint bash  $job_script \n";
        }
        else {
                print STCH " bash  $job_script \n";
        }

        close STCH;

    $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 4a: Perform assembly using metaSPAdes: both paired end and single end reads in the same step##########
sub metaSPAdes_assembly{
        $current_job_file1 = "j4_".$sample_name."_metaSPAdes_assembly.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage4\n";
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
        
	if (!$step_number or $step_number == -1) {
                chomp $last_jobid1;
                print STCH "#SBATCH --depend=afterok:$last_jobid1 \n";
        }
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  

        print STCH "OUTDIR=".$assembly_dir_metaSPAdes."\n";
        print STCH "IN1=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n\n";
        print STCH "IN2=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n\n";
        print STCH "INs=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n\n";

        print STCH "if [ -e \${IN} ] \n";
        print STCH "then\n";
        print STCH "    if [ ! -e $status_log/j4_finished_metaSPAdes_assembly\${SLURM_ARRAY_TASK_ID} ]\n";
        print STCH "    then\n";
        
	# If assembly directory already exists, the assembly failed previously. Need to remove old assembly directory.
        print STCH "            if [ -d \${OUTDIR} ]\n";
        print STCH "            then\n";
        print STCH "                    rm -rf \${OUTDIR} \n";
        print STCH "            fi\n";
        print STCH "            $metaSPAdes --memory $SPAdes_mem -t \$SLURM_CPUS_PER_TASK -1 \${IN1} -2 \${IN2} -o \${OUTDIR} \n";
        print STCH "            OUT=\$?\n";
        print STCH '            if [ ${OUT} -ne 0  ]',"\n"; # did not finish successfully
        print STCH "            then\n";
        print STCH "                    echo \"Fatal Error : MetaSPAdes did not finish correctly.\"  \n";
        print STCH "                    exit 1\n";
        print STCH "            fi\n";
        print STCH "            cp \${OUTDIR}/contigs.fasta $seq_dir/spades_contigs.fasta\n";
        print STCH "            touch $status_log/j4_finished_metaSPAdes_assembly\${SLURM_ARRAY_TASK_ID}  \n";
        print STCH "    fi\n";
        print STCH "fi\n";
        close STCH;

        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 4b: Perform long-read OR hybrid assembly using dragonflye
sub dragonflye_assembly{
        $current_job_file1 = "j4b_".$sample_name."_dragonflye_assembly.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage4b\n";
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
        
	if (!$step_number or $step_number == -1) {
                chomp $last_jobid1;
                print STCH "#SBATCH --depend=afterok:$last_jobid1 \n";
        }
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate dragonflye\n";  

        print STCH "OUTDIR=".$assembly_dir_dragonflye."\n";
        print STCH "IN1=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n\n";
        print STCH "IN2=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n\n";
        print STCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.LR.fastq\n";
	
	my $dragon_com1="dragonflye --force --tmpdir /dev/shm --racon 4 --outdir \$OUTDIR --reads \$IN\n";
	my $dragon_com2="dragonflye --force --tmpdir /dev/shm --racon 4 --outdir \$OUTDIR --reads \$IN --R1 \$IN1 --R2 \$IN2\n";
	my $dragon_com="";
	if ($assembly_type eq "h"){
		$dragon_com = $dragon_com2;
	}
	else{
		$dragon_com = $dragon_com1;
	}

        print STCH "if [ -e \${IN} ] \n";
        print STCH "then\n";
        print STCH "    if [ ! -e $status_log/j4b_finished_dragonflye_assembly\${SLURM_ARRAY_TASK_ID} ]\n";
        print STCH "    then\n";
        
	# If assembly directory already exists, the assembly failed previously. Need to remove old assembly directory.
        print STCH "            if [ -d \${OUTDIR} ]\n";
        print STCH "            then\n";
        print STCH "                    rm -rf \${OUTDIR} \n";
        print STCH "            fi\n";

        print STCH "            $dragon_com";
        print STCH "            OUT=\$?\n";
        print STCH '            if [ ${OUT} -ne 0  ]',"\n"; # did not finish successfully
        print STCH "            then\n";
        print STCH "                    echo \"Fatal Error : dragonflye did not finish correctly.\"  \n";
        print STCH "                    exit 1\n";
        print STCH "            fi\n";
        print STCH "            cp \${OUTDIR}/flye.fasta $seq_dir/dragonflye_contigs.fasta\n";
        print STCH "            touch $status_log/j4b_finished_dragonflye_assembly\${SLURM_ARRAY_TASK_ID}  \n";
        print STCH "    fi\n";
        print STCH "fi\n";
        close STCH;

        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 4c: Perform hybrid assembly using unicycler
sub unicycler_assembly{
        $current_job_file1 = "j4c_".$sample_name."_unicycler_assembly.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage4c\n";
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
        
	if (!$step_number or $step_number == -1) {
            chomp $last_jobid1;
             print STCH "#SBATCH --depend=afterok:$last_jobid1 \n";
        }
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  

        print STCH "OUTDIR=".$assembly_dir_unicycler."\n";
        print STCH "IN1=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n\n";
        print STCH "IN2=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n\n";
        print STCH "INL=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.LR.fastq\n";

        print STCH "if [ -e \${IN} ] \n";
        print STCH "then\n";
        print STCH "    if [ ! -e $status_log/j4c_finished_unicycler_assembly\${SLURM_ARRAY_TASK_ID} ]\n";
        print STCH "    then\n";
        
	# If assembly directory already exists, the assembly failed previously. Need to remove old assembly directory.
        print STCH "            if [ -d \${OUTDIR} ]\n";
        print STCH "            then\n";
        print STCH "                    rm -rf \${OUTDIR} \n";
        print STCH "            fi\n";

        print STCH "            unicycler -1 \$IN1 -2 \$IN2 -l \$INL --min_fasta_length 1000 -t $num_CPU --mode normal --out \$OUTDIR \n";
        print STCH "            OUT=\$?\n";
        print STCH '            if [ ${OUT} -ne 0  ]',"\n"; # did not finish successfully
        print STCH "            then\n";
        print STCH "                    echo \"Fatal Error : unicycler did not finish correctly.\"  \n";
        print STCH "                    exit 1\n";
        print STCH "            fi\n";
        print STCH "            cp \${OUTDIR}/assembly.fasta $seq_dir/unicycler_contigs.fasta\n";
        print STCH "            touch $status_log/j4c_finished_unicycler_assembly\${SLURM_ARRAY_TASK_ID}  \n";
        print STCH "    fi\n";
        print STCH "fi\n";
        close STCH;

        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 4d: Perform assembly using SPAdes: both paired end and single end reads in the same step##########
sub SPAdes_assembly{
        $current_job_file1 = "j4d_".$sample_name."_SPAdes_assembly.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
        print STCH "#SBATCH --job-name=${sample_path}:stage4\n";
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";

        if (!$step_number or $step_number == -1) {
                chomp $last_jobid1;
                print STCH "#SBATCH --depend=afterok:$last_jobid1 \n";
        }
		print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";

        print STCH "OUTDIR=".$assembly_dir_metaSPAdes."\n";
        print STCH "IN1=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n\n";
        print STCH "IN2=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n\n";
        print STCH "INs=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n\n";

        print STCH "if [ -e \${IN} ] \n";
        print STCH "then\n";
        print STCH "    if [ ! -e $status_log/j4d_finished_SPAdes_assembly\${SLURM_ARRAY_TASK_ID} ]\n";
        print STCH "    then\n";

        # If assembly directory already exists, the assembly failed previously. Need to remove old assembly directory.
        print STCH "            if [ -d \${OUTDIR} ]\n";
        print STCH "            then\n";
        print STCH "                    rm -rf \${OUTDIR} \n";
        print STCH "            fi\n";
        #TODO: check for empty singleton file to avoid SPAdes crash
        print STCH "            $SPAdes -t \$SLURM_CPUS_PER_TASK --memory $SPAdes_mem -1 \${IN1} -2 \${IN2} -o \${OUTDIR} \n";
        print STCH "            OUT=\$?\n";
        print STCH '            if [ ${OUT} -ne 0  ]',"\n"; # did not finish successfully
        print STCH "            then\n";
        print STCH "                    echo \"Fatal Error : SPAdes did not finish correctly.\"  \n";
        print STCH "                    exit 1\n";
        print STCH "            fi\n";
        print STCH "            cp \${OUTDIR}/contigs.fasta $seq_dir/spades_contigs.fasta\n";
        print STCH "            touch $status_log/j4d_finished_SPAdes_assembly\${SLURM_ARRAY_TASK_ID}  \n";
        print STCH "    fi\n";
        print STCH "fi\n";
        close STCH;

        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 5a: Map reads back to contigs for short reads##########
sub map_to_assembly{
        # this is the job script
        my $job_script = $job_files_dir."/j5_".$sample_name."_map_to_assembly_script.sh";
        open(JSTCH, ">$job_script") or die $!;
        chmod 0755, $job_script;
        print JSTCH "#!/bin/bash\n";
        print JSTCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print JSTCH "conda activate vs\n";  
        print JSTCH "IN1=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n";
        print JSTCH "IN2=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n";
        print JSTCH "INs=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n";
        #print JSTCH "REF=".$assembly_dir_metaSPAdes."/contigs.fasta\n";
        print JSTCH "OUTFILE1=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.mapped_pe.fastq.gz\n";
        print JSTCH "OUTFILE2=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.fastq.gz\n";
        print JSTCH "OUTFILE3=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.mapped_se.fastq.gz\n";
        print JSTCH "OUTFILE4=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_se.fastq.gz\n";
        print JSTCH "OUTFILE5=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.R1.fastq.gz\n";
        print JSTCH "OUTFILE6=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.R2.fastq.gz\n";
        print JSTCH "OUTFILE7=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.se.fastq.gz\n";
        print JSTCH "TIMEFILE=$sample_dir"."/".$sample_name.".QC.assembly.map_time.txt\n";
		#set REF depending on which type of assembly you are using
		if ($assembly_type eq "s"){
		print JSTCH "REF=".$assembly_dir_metaSPAdes."/contigs.fasta\n";
		}
		if ($assembly_type eq "h" and $assembly_mode eq "i"){
		print JSTCH "REF=".$assembly_dir_unicycler."/assembly.fasta\n";
		}
		if ($assembly_type eq "h" and $assembly_mode eq "m"){
		print JSTCH "REF=".$assembly_dir_dragonflye."/flye.fasta\n";
		}
        # Check if mapping finished successfully. If true, resubmitting will skip this step and keep results
        print JSTCH "if [ ! -e $status_log/j5_map_to_assembly_finished ]\n";
        print JSTCH "then\n";
        print JSTCH "   date > \${TIMEFILE}\n";
        print JSTCH "   bbmap.sh in1=\${IN1} in2=\${IN2} ref=\${REF} covstats=$sample_dir"."/bbmap_assembly_PE_covstats.txt statsfile=$sample_dir"."/bbmap_assembly_PE_statsfile.txt fast=t nodisk=t ow 32bit=t pigz=t -Xmx${slurm_mem} outm=\${OUTFILE1} outu=\${OUTFILE2} \n";
        print JSTCH "   bbmap.sh in=\${INs}  ref=\${REF} covstats=$sample_dir"."/bbmap_assembly_SE_covstats.txt statsfile=$sample_dir"."/bbmap_assembly_SE_statsfile.txt fast=t nodisk=t ow 32bit=t pigz=t -Xmx${slurm_mem} outm=\${OUTFILE3} outu=\${OUTFILE4} \n";
        print JSTCH "   $reformat in=\${OUTFILE2} out1=\${OUTFILE5} out2=\${OUTFILE6} outs=\${OUTFILE7} ow\n";
        print JSTCH "   OUT=\$?\n";
        print JSTCH '   if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
        print JSTCH "   then\n";
        print JSTCH "           echo \"Fatal Error trying to run bbmap.\"  \n";
        print JSTCH "           exit 1\n";
        print JSTCH "   else\n";
        print JSTCH "           touch $status_log/j5_map_to_assembly_finished\n";
        print JSTCH "   fi\n";

        #if mapping finished succussfully, make a file named finished
        print JSTCH "   date >> \${TIMEFILE}\n";
        print JSTCH "fi\n\n";
        close JSTCH;

        # here to submit
        $current_job_file1 = "j5_".$sample_name."_map_to_assembly.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage5\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n"; # request total memory of 128G
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
        if (!$step_number or $step_number == -1) {
                print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        }
        if ($use_checkpoint) {
                print STCH "checkpoint bash  $job_script \n";
        }
        else {
                print STCH " bash  $job_script \n";
        }

        close STCH;

    $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 5b: Map reads back to contigs for long reads##########
sub map_to_assembly_minimap{
        # this is the job script
        my $job_script = $job_files_dir."/j5b_".$sample_name."_minimap_to_assembly_script.sh";
        open(JSTCH, ">$job_script") or die $!;
        chmod 0755, $job_script;
        print JSTCH "#!/bin/bash\n";
        print JSTCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.LR.fastq\n";
        print JSTCH "OUTFILE1=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.assembly.minimap.sam\n";
        print JSTCH "OUTFILE11=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.assembly.minimap_sorted.bam\n";
        print JSTCH "OUTFILE2=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.assembly.mapped.LR.fastq\n";
        print JSTCH "OUTFILE3=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.assembly.unmapped.LR.fastq\n";
        print JSTCH "TIMEFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.mapL_time.txt\n";
        print JSTCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print JSTCH "conda activate vs\n";

		#set REF depending on which type of assembly you are using
		if ($assembly_type eq "l"){
			print JSTCH "REF=".$assembly_dir_dragonflye."/flye.fasta\n";
		}
		if ($assembly_type eq "h" and $assembly_mode eq "i"){
			print JSTCH "REF=".$assembly_dir_unicycler."/assembly.fasta\n";
		}
                if ($assembly_type eq "h" and $assembly_mode eq "m"){
                	print JSTCH "REF=".$assembly_dir_dragonflye."/flye.fasta\n";
                }
        # Check if mapping finished successfully. If true, resubmitting will skip this step and keep results
        print JSTCH "if [ ! -e $status_log/j5b_map_to_assembly_finished ]\n";
        print JSTCH "then\n";
        print JSTCH "   date > \${TIMEFILE}\n";
        print JSTCH "   minimap2 -ax map-ont \$REF \$IN > \$OUTFILE1\n";
        print JSTCH "   OUT=\$?\n";
        print JSTCH '   if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
        print JSTCH "   then\n";
        print JSTCH "           echo \"Fatal Error trying to run minimap2.\"  \n";
        print JSTCH "           exit 1\n";
        print JSTCH "   else\n";
        print JSTCH "		samtools fastq -F 4  \$OUTFILE1 > \$OUTFILE2\n";
		print JSTCH "		samtools fastq -f 4 \$OUTFILE1 > \$OUTFILE3\n";
        print JSTCH "           samtools sort \$OUTFILE1 -O bam -o \$OUTFILE11\n";
        print JSTCH "           samtools index \$OUTFILE11\n";
        print JSTCH "           samtools idxstats \$OUTFILE11 > minimap__assembly_LR_covstats.txt\n";  
        print JSTCH "           touch $status_log/j5b_map_to_assembly_finished\n";
        print JSTCH "   fi\n";

        #if mapping finished succussfully, make a file named finished
        print JSTCH "   date >> \${TIMEFILE}\n";
        print JSTCH "fi\n\n";
        close JSTCH;

        # here to submit
        $current_job_file1 = "j5b_".$sample_name."_map_to_assembly.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
		print STCH "#SBATCH --job-name=${sample_path}:stage5b\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n"; 
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
        if (!$step_number or $step_number == -1) {
                print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        }
        if ($use_checkpoint) {
                print STCH "checkpoint bash  $job_script \n";
        }
        else {
                print STCH " bash  $job_script \n";
        }

        close STCH;

    $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}


##########STAGE 4b/5c: Skip assembly/contig mapping##########
sub skip_assembly{
        $current_job_file1 = "j4_".$sample_name."_skip_assembly.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage4\n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";

        if (!$step_number or $step_number == -1) {
                print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        }
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  

        print STCH "IN1=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n\n";
        print STCH "IN2=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n\n";
        print STCH "INs=".$sample_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n\n";
        print STCH "OUTsingle1=".$seq_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n\n";
        print STCH "OUTsingle2=".$seq_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.se.fastq.gz\n\n";
        print STCH "OUTsingle3=".$seq_dir."/${sample_name}.RemoveAdapter.bbduk.RefGenome.unmapped.se.fasta\n\n";

	# Check whether skip assembly finished successfully. If true, resubmitting will skip this step and keep results
        print STCH "OUT=1\n";
        print STCH "if [ ! -e $status_log/j4_skip_assembly_finished ]\n";
        print STCH "then\n";
        print STCH "    $seqtk seq -A \${INs} > \${OUTsingle3}\n";
        print STCH "    OUT=\$?\n";
        print STCH '    if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
        print STCH "    then\n";
        print STCH "            echo \"Fatal Error trying to assembly.\"  \n";
        print STCH "            exit 1\n";

        # if finished, make a file named _finished
        print STCH "    else\n";
        print STCH "            touch $status_log/j4_skip_assembly_finished\n";
        print STCH "    fi\n";
        print STCH "fi\n\n";
        close STCH;

        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 6: Stitch paired reads together if unmapped to assembly##########
sub stitching{
	$current_job_file1 = "j6_".$sample_name."_stitching.sh";
	$current_job_file2 = $current_job_file1;
		
	######################################################################
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage6\n";
	print STCH "#SBATCH --mem=${slurm_mem}\n";
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".err\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  
        if ($VStrack eq "V") {
            print STCH "IN1=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe1.fastq.gz\n";
            print STCH "IN2=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe2.fastq.gz\n";
            print STCH "OUTFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe\n";
            #my $full_sample_name=$sample_dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe\n";
        }else {
            print STCH "IN1=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.R1.fastq.gz\n";
            print STCH "IN2=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.R2.fastq.gz\n";
            print STCH "OUTFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.pear\n";
            #my $full_sample_name=$sample_dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe\n";
        }
        print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".stitching_time.txt\n";
	print STCH "OVERLAPLENGTH=$sample_dir"."/".$sample_name.".stitching_OVERLAPLENGTH.txt\n";
	my $stitching_report_file = $sample_dir."/".$sample_name.".stitching.report";

	#Check whether stitching finished successfully. If true, resubmitting will skip this step and keep results.
	print STCH "if [ ! -e $status_log/j6_stitching_finished ]\n";
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n\n";
	print STCH "	$pear -j \${SLURM_CPUS_PER_TASK} -v 10 -f \${IN1} -r \${IN2} -o \${OUTFILE} > $stitching_report_file\n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to run pear.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n\n";

	# Pear produces 4 files: assembled, discarded, unassembled.forward, and unassembled.reverse. This step joins all non-discarded reads into a single file
        print STCH "	".$run_script_path."Combine_joined_unjoin_reads.pl \${OUTFILE}\n";
        print STCH "	OUT2=\$?\n";
        print STCH '	if [ ${OUT2} -ne 0 ]',"\n"; # did not finish successfully
        print STCH "	then\n";
        print STCH "		echo \"Fatal Error trying to post process pear files.\"  \n";
        print STCH "		exit 1\n"; 
        print STCH "	fi\n\n";

	# Compress fastq files
        print STCH "	$pigz --force \${OUTFILE}.*\n";
        print STCH "	OUT3=\$?\n";
        print STCH '	if [ ${OUT3} -ne 0 ]',"\n"; # did not finish successfully
        print STCH "	then\n";
        print STCH "		echo \"Fatal Error trying to compress pear files.\"  \n";
        print STCH "		exit 1\n"; 
        print STCH "	fi\n\n";
	
        # Stitching finished successfully, make the log file 
	print STCH "	touch $status_log/j6_stitching_finished\n";
	print STCH "	date >> \${TIMEFILE}\n";
	print STCH "fi\n\n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 7: Combine sequences into single file##########
sub get_all_seqs{
        $current_job_file1 = "j7_".$sample_name."_get_all_seqs.sh";
        $current_job_file2 = $current_job_file1;

        ######################################################################
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage7\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  
	
        if ($VStrack eq "D") {
        print STCH "INsingle1=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.se.fastq.gz\n";
		print STCH "INsingle2=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_se.fastq.gz\n";
        print STCH "OUTsingle1=$seq_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.se.fasta\n";
        print STCH "OUTsingle2=$seq_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_se.fasta\n";
        print STCH "INstitch=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.pear.stitched.fastq.gz\n";
        print STCH "OUTstitch=$seq_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.stitched.fasta\n";
        
        print STCH "INlong=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.assembly.unmapped.LR.fastq\n";
        print STCH "OUTlong=$seq_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.assembly.unmapped.LR.fasta\n";

        }

        else {
		print STCH "INstitch=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe.pear.stitched.fastq.gz\n";
		print STCH "OUTstitch=$seq_dir"."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.pe.stitched.fasta\n";
		print STCH "INlong=$sample_dir"."/".$sample_name.".RemoveAdapter.fastp.RefGenome.unmapped.LR.fastq.gz\n";
        }

        print STCH "OUTFILE=$sample_dir/${sample_name}.allSeqs.fasta\n\n";
		print STCH "if [ ! -e $status_log/j7_get_all_seqs_finished ]\n";
		print STCH "then\n";


	# convert fastq to fasta
		if ($assembly_type eq "s" or $assembly_type eq "h"){
			if ($VStrack eq "D") {
			print STCH "    $seqtk seq -A \${INsingle1} > \${OUTsingle1}\n";
			print STCH "    $seqtk seq -A \${INsingle2} > \${OUTsingle2}\n";
			}
			print STCH "    $seqtk seq -A \${INstitch} > \${OUTstitch}\n";
		}
		if ($assembly_type eq "l" or $assembly_type eq "h"){
			print STCH "    $seqtk seq -A \${INlong} > \${OUTlong}\n";
		}
	
	# combine files 	
	print STCH "	cat $seq_dir/*.fasta >> \${OUTFILE} \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0  ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error in combining sequences.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	touch $status_log/j7_get_all_seqs_finished \n";
	print STCH "fi\n";
	close STCH;
    
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 8: This step will cluster sequences using mmseqs to reduce redundancy##########
sub run_mmseqs{
	$current_job_file1 = "j8_".$sample_name."_cluster.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage8\n";
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --mem=${slurm_mem}\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
        
	print STCH "if [ ! -e ${status_log}/j8_cluster_finished ]\n";
	print STCH "then\n";
	print STCH "	CLUSTER_IN=${sample_dir}/${sample_name}.allSeqs.fasta\n";
	print STCH "	CLUSTER_OUT=${sample_dir}"."/".${sample_name}.".allSeqs.mmseqs\n";

	print STCH "	if [ ! -e ${status_log}/j8_cluster_finished ]\n"; 
	print STCH "	then\n";
	print STCH " 		$mmseqs_easylinclust \${CLUSTER_IN} \${CLUSTER_OUT} $mmseqs_tmpdir $mmseqs_parameter\n";
	print STCH "		OUT=\$?\n"; 
	print STCH "		if [ \${OUT} -ne 0 ]\n"; # did not finish successfully
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to run mmseqs.\"  \n";
	print STCH "			exit 1\n"; 
	# if finished successfully, make a file named cluster_finished
	print STCH "		else\n";
	print STCH "			touch $status_log/j8_cluster_finished\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;

    $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 9: Split unique reads to small files for RepeatMasker##########
sub pre_RepeatMasker{
	$current_job_file1 = "j9_".$sample_name."_pre_RepeatMasker.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage9\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "IN=".${sample_name}.".allSeqs.mmseqs_rep_seq.fasta\n";
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "RM_DIR=".$repeatmasker_dir."\n\n";

	print STCH "if [ ! -e $status_log/j9_pre_RepeatMasker_finished ]\n";
	print STCH "then\n";

	# if the directory already exist, remove it.
	print STCH "	if [ -d \${RM_DIR} ]\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${RM_DIR} \n";
	print STCH "	fi\n";

	# make repeatmasker directory
	print STCH "	mkdir -p \${RM_DIR} \n\n";

	# split into smaller files
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i  \${IN} -o $repeatmasker_dir -f $file_number_of_RepeatMasker -n $num_seq_per_file_Repeatmasker \n";
	print STCH "	CHECK1=\$?\n"; # if finishes value is 0
	print STCH "	".$run_script_path."check_split_RepeatMasker.pl \${SAMPLE_DIR}\n";
	print STCH "	CHECK2=\$?\n"; # if finishes correctly value is 0
	print STCH "	if [ \${CHECK1} -ne 0 -o \${CHECK2} -ne 0 ] \n"; 
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error at preparing for RepeatMasker.\"  \n";
	print STCH "		exit 1\n"; 
	# if finished succussfully, make log file
	print STCH "	else\n";
	print STCH "		touch $status_log/j9_pre_RepeatMasker_finished \n";
	print STCH "	fi\n";
	print STCH "fi\n\n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 10: Run RepeatMasker###########
sub RepeatMasker{
	# this is the job script
	my $job_script = $job_files_dir."/j10_".$sample_name."_RepeatMasker_script.sh";
	open(JSTCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print JSTCH "#!/bin/bash\n";         
	print JSTCH "RMOUT=$repeatmasker_dir/${sample_name}.allSeqs.mmseqs_rep_seq.fasta_file\${SLURM_ARRAY_TASK_ID}.fasta.masked\n";#
	print JSTCH "RMIN=$repeatmasker_dir/${sample_name}.allSeqs.mmseqs_rep_seq.fasta_file\${SLURM_ARRAY_TASK_ID}.fasta\n";#
	print JSTCH "RMOTHER=$repeatmasker_dir/${sample_name}.allSeqs.mmseqs_rep_seq.fasta_file\${SLURM_ARRAY_TASK_ID}.fasta.out\n";
	print JSTCH "RM_dir=".$repeatmasker_dir."\n\n";
	print JSTCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print JSTCH "conda activate vs\n";  
	
	print JSTCH 'if [ -e $RMIN ]',"\n"; # input file exist
	print JSTCH "then\n";
	print JSTCH '	if [ ! -e $RMOTHER ]',"\n"; # don't have RepeatMasker output ".out" file, means RepeatMasker never ran or finished
	print JSTCH "	then\n";
	print JSTCH "		$repeat_masker -qq -par \$SLURM_CPUS_PER_TASK \$RMIN \n"; # run RepeatMasker, -par number of processors to use
	print JSTCH '     	CHECK=$?',"\n"; # if finished correctly the value will be 0
	print JSTCH "		if [ \${CHECK} -ne 0 ] \n"; 
	print JSTCH "		then\n";
	print JSTCH "			echo \"Fatal Error at RepeatMasker.\"  \n";
	print JSTCH "			exit 1\n"; 
	# if finished succussfully, make log file
	print JSTCH "		fi\n";
	print JSTCH "	fi\n\n";

 	print JSTCH '	if [ ! -e $RMOUT ]',"\n"; #sometimes repeatmasker does not find any repeat in input files, in these cases no .masked file will be generated.
	print JSTCH "	then\n";
	print JSTCH '		cp ${RMIN} ${RMOUT}',"\n";
	print JSTCH "	fi\n";
	print JSTCH "fi\n";
	close JSTCH;

	# here to submit
	$current_job_file1 = "j10_".$sample_name."_RepeatMasker.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage10\n";
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --array=1-$file_number_of_RepeatMasker\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	#print STCH "set -x\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 11: Check RepeatMasker output##########
sub post_RepeatMasker{
	$current_job_file1 = "j11_".$sample_name."_post_RepeatMasker.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage11\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".RepeatMasker_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j19_post_RepeatMasker_finished ]\n";
	print STCH "then\n";
	print STCH "	".$run_script_path."RepeatMasker_SequenceQualityControl.pl ".$sample_dir."\n";
	print STCH '	CHECK1=$?',"\n"; # if finishes value is 0
	print STCH "	".$run_script_path."RepeatMasker_check_SequenceQualityControl.pl \${SAMPLE_DIR} > \${QC_Record} \n";
	print STCH '	CHECK2=$?',"\n"; # if finishes correctly value is 0
	print STCH '	if [ ${CHECK1} -ne 0 -o ${CHECK2} -ne 0 ]',"\n"; # either one is not finished correctly
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error at RepeatMasker SequenceQualityControl.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	else\n";
	print STCH "		touch $status_log/j19_post_RepeatMasker_finished \n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 12: Split input file for MMseqs against viral-only database##########
sub split_for_MMseqs_VIRUSDB{
	$current_job_file1 = "j12_".$sample_name."_MMseqs_VIRUSDB_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage12\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "MMseqs_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$MMseqs_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "IN=".$sample_name.".goodSeq.fa\n\n";

	print STCH "if [ ! -e $status_log/j12_MMseqs_VIRUSDB_split_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	if [ -d \${MMseqs_VIRUSDB_DIR} ] \n"; # If the dir already exists, remove.
	print STCH "	then\n";
	print STCH "		rm -rf \${MMseqs_VIRUSDB_DIR}\n";
	print STCH "	fi\n\n";

	# make the directory
	print STCH "	mkdir -p \${MMseqs_VIRUSDB_DIR}\n";
	
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${IN} -o \${MMseqs_VIRUSDB_DIR} -f 1 -n $num_seq_per_file_MMseqs_VIRUSDB \n";
	print STCH '	CHECK=$?',"\n";
	print STCH '	if [ ${CHECK} -ne 0 ]',"\n"; 
	print STCH "	then\n";
	# split to -n number of files, this number should be consistent with the number of 
	# MMseqs job array submitted bellow
	print STCH "		echo \"Fatal Error trying to split for MMseqs_VIRUSDB\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	else\n";
	print STCH "		touch $status_log/j12_MMseqs_VIRUSDB_split_finished \n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;

	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 13: Run mmseqs2 against viral-only database##########
sub submit_job_array_MMseqs_VIRUSDB{
        # this is the job script
        my $job_script = $job_files_dir."/j13_".$sample_name."_mmseqs_VIRUSDB_run_script.sh";
        open(STCH, ">$job_script") or die $!;
        chmod 0755, $job_script;

        print STCH "#!/bin/bash\n";         
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";
        print STCH "MMseqs_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$MMseqs_VIRUSDB_DIR_SUFFIX."\n";
        print STCH "MMseqsOUT=",'${MMseqs_VIRUSDB_DIR}',"/",$sample_name.".goodSeq.fa_file".'${SLURM_ARRAY_TASK_ID}',".mmseqsv.out\n"; # full path
        print STCH "QUERY=",'${MMseqs_VIRUSDB_DIR}',"/".$sample_name.".goodSeq.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n";
        print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".mmseqs_virusDB_time.txt\n";

        print STCH 'if [ -s $QUERY ]',"\n"; #check if the file is empty
        print STCH "then\n";
        #if the output file does not exist, run and check the completeness of the output file
        print STCH '    if [ ! -e $MMseqsOUT ]',"\n";
        print STCH "    then\n";
        print STCH "		date > \${TIMEFILE}\n";
        print STCH "            $mmseqs_easysearch \${QUERY} $db_MMSEQS_NT_VIRUS \${MMseqsOUT} $mmseqs_tmpdir --max-accept 25 --threads \$SLURM_CPUS_PER_TASK --start-sens 1 --sens-steps 3 -s 7 -e 1.0E-4 --split-memory-limit \${SLURM_MEM_PER_NODE} --search-type 3 --format-output \"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen\"","\n";
        print STCH '            CHECK=$?',"\n";
        print STCH '      	if [ ${CHECK} -ne 0 ]',"\n";
        print STCH "            then\n";
        print STCH "                    echo \"Fatal Error : mmseqs easy-search did not finish correctly.\"  \n";
        print STCH "                    exit 1\n";
        print STCH "            else\n";
#        print STCH "               echo 'MMSEQS_EASYSEARCH_COMPLETED' >> \${MMseqsOUT}","\n";
        print STCH "            	date >> \${TIMEFILE}\n";
      	print STCH "			touch $status_log/j13_finished_mmseqs_easysearch \n";
        print STCH "            fi\n";

        #if the output file exists, check the completeness of the output file
        print STCH "    else\n";

        print STCH "		if [ ! -e $status_log/j13_finished_mmseqs_easysearch ]\n";
        print STCH "            then\n";
        print STCH "            	date > \${TIMEFILE}\n";
        print STCH "                    $mmseqs_easysearch \${QUERY} $db_MMSEQS_NT_VIRUS \${MMseqsOUT} $mmseqs_tmpdir --max-accept 25 --threads \$SLURM_CPUS_PER_TASK --start-sens 1 --sens-steps 3 -s 7 -e 1.0E-4 --split-memory-limit \${SLURM_MEM_PER_NODE} --search-type 3 --format-output \"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen\"","\n";
        print STCH '                    CHECK=$?',"\n";
        print STCH '                    if [ ${CHECK} -ne 0 ]',"\n";
        print STCH "                    then\n";
        print STCH "                            echo \"Fatal Error : mmseqs easy-search did not finish correctly.\"  \n";
        print STCH "                            exit 1\n";
        print STCH "                    else\n";
        print STCH "                    	date >> \${TIMEFILE}\n";
    	print STCH "				touch $status_log/j13_finished_mmseqs_easysearch \n";   
        print STCH "                    fi\n";
        print STCH "            fi\n";
        print STCH "    fi\n";
        print STCH "fi";
        close STCH;

        # here to submit
        $current_job_file1 = "j13_".$sample_name."_MMseqs_VIRUSDB_run.sh";
        open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
        print STCH "#!/bin/bash\n";  
        print STCH "#SBATCH --partition=$partition\n";
        print STCH "#SBATCH --account=$account\n";
        print STCH "#SBATCH --job-name=${sample_path}:stage13\n";
        print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
        print STCH "#SBATCH --mem=${slurm_mem}\n";
        print STCH "#SBATCH --array=1-$file_number_of_MMseqs_VIRUSDB\n";
        if (!$step_number or $step_number == -1) {
                print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
        }
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
        print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
        if ($use_checkpoint) {
                print STCH "checkpoint bash  $job_script \n";
        }
	else {
              	print STCH "bash $job_script \n";
        }

	close STCH;
        $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 14: Parse output of mmseqs against viral-only database##########
sub parse_MMseqs_VIRUSDB {
	#my $parser_path = "/export/virusseeker/scripts/qkt_scripts/blast_change/BLAST_VIRUSDB_parser.py";
	$current_job_file1 = "j14_".$sample_name."_Parse_MMseqs_VIRUSDB.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n"; 
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage14\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "#SBATCH --array=1-$file_number_of_MMseqs_VIRUSDB\n";
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "MMseqs_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$MMseqs_VIRUSDB_DIR_SUFFIX."/\n";
	print STCH "MMseqsOUT=",$sample_name.".goodSeq.fa_file".'${SLURM_ARRAY_TASK_ID}',".mmseqsv.out\n";#name only, not full path
	print STCH "MMseqsIN=",$sample_name.".goodSeq.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n";#full path
	print STCH "PARSED=",'${MMseqs_VIRUSDB_DIR}',"/".$sample_name.".goodSeq.fa_file".'${SLURM_ARRAY_TASK_ID}',".mmseqsv.parsed\n";
	print STCH "MMseqs_VIRUSDB_Filtered_fa=",'${MMseqs_VIRUSDB_DIR}',"/".$sample_name.".goodSeq.fa_file".'${SLURM_ARRAY_TASK_ID}',".MMseqs_VIRUSDB_filtered.fa\n\n";

	print STCH 'if [ -f ${MMseqs_VIRUSDB_DIR}/${MMseqsOUT} ]',"\n"; 
	print STCH "then\n";

	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -f $PARSED ]',"\n";
	print STCH "	then\n";
	print STCH "            python ".$bash_bin_path."/parse_blast_Virus_noBP.py -d \${MMseqs_VIRUSDB_DIR}/ -i \${MMseqsOUT} -f \${MMseqsIN}\n";
	print STCH "		".$run_script_path."check_MMseqs_parsed_file.pl \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed, if not completed, exit
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${MMseqsOUT}.\"  \n";
	print STCH '			exit 1',"\n";
	print STCH "		fi\n";

	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_MMseqs_parsed_file.pl  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not correctly completed exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "		then\n"; # rerun the parser
	print STCH "			python ".$bash_bin_path."/parse_blast_Virus_noBP.py -d \${MMseqs_VIRUSDB_DIR}/ -i \${MMseqsOUT} -f \${MMseqsIN}\n";
	print STCH "			".$run_script_path."check_MMseqs_parsed_file.pl  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed, if not completed, exit
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${MMseqsOUT}.\"  \n";
	print STCH '				exit 1',"\n";
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 15: Split input file for blastx against viral-only database##########
sub split_for_BLASTX_VIRUSDB {
	$current_job_file1 = "j15_".$sample_name."_BLASTX_VIRUSDB_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage15\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  

	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "BLASTX_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "MMseqs_VIRUSDB_Filtered_fa=".$sample_name.".MMseqs_VIRUSDB_Filtered.fa\n";
	print STCH "MMseqs_VIRUSDB_Hit_fa=".$sample_name.".MMseqs_VIRUSDB_HIT.fa\n";
	print STCH "MMseqs_VIRUSDB_Phage_fa=".$sample_name.".MMseqs_VIRUSDB_Phage.fa\n";
	print STCH "MMseqs_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$MMseqs_VIRUSDB_DIR_SUFFIX."\n\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".MMseqs_VIRUSDB_finish_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j15_finished_split_for_BLASTX_VIRUSDB ]\n";
	print STCH "then\n";
	# check to make sure the number of reads add togetehr equals the total input reads
	print STCH "	".$run_script_path."check_MMseqs_VIRUSDB_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : MMseqs_VIRUSDB did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	# if the BLASTX Virus database directory already exist, remove it
	print STCH '	if [ -d $BLASTX_VIRUSDB_DIR ] ',"\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${BLASTX_VIRUSDB_DIR}\n";
	print STCH "	fi\n";

	# make BLASTX Virus database directory
	print STCH "	mkdir -p \${BLASTX_VIRUSDB_DIR}\n";

	# pool all the MMseqs VIRUSDB hits into a single file
	print STCH "	if [ -e $sample_dir/\${MMseqs_VIRUSDB_Hit_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${MMseqs_VIRUSDB_Hit_fa}\n";
	print STCH "	fi\n";
	print STCH "	cat \${MMseqs_VIRUSDB_DIR}/*.MMseqs_VIRUSDB_hit.fa >> $sample_dir/\${MMseqs_VIRUSDB_Hit_fa}\n\n";

	# pool all the MMseqs VIRUSDB filtered reads into a single file
	print STCH "	if [ -e $sample_dir/\${MMseqs_VIRUSDB_Filtered_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${MMseqs_VIRUSDB_Filtered_fa}\n";
	print STCH "	fi\n";
	print STCH "	cat \${MMseqs_VIRUSDB_DIR}/*.MMseqs_VIRUSDB_filtered.fa >> $sample_dir/\${MMseqs_VIRUSDB_Filtered_fa}\n\n";

	# pool all the MMseqs VIRUSDB phage hits into a single file
	print STCH "	if [ -e $sample_dir/\${MMseqs_VIRUSDB_Phage_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${MMseqs_VIRUSDB_Phage_fa}\n";
	print STCH "	fi\n";
	print STCH "	cat \${MMseqs_VIRUSDB_DIR}/*.MMseqs_VIRUSDB_phage.fa >> $sample_dir/\${MMseqs_VIRUSDB_Phage_fa}\n\n";

	# split into smaller files
	# split to -n number of files, this number should be consistent with 
	# the number of diamond job array submitted bellow
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d \${SAMPLE_DIR} -i  \${MMseqs_VIRUSDB_Filtered_fa} -o \${BLASTX_VIRUSDB_DIR} -f $file_number_of_BLASTX_VIRUSDB  -n $num_seq_per_file_BLASTX_VIRUSDB \n";
	print STCH '	CHECK1=$?',"\n";
	print STCH "	".$run_script_path."check_split_BLASTX_VIRUSDB.pl \${SAMPLE_DIR}\n";
	print STCH '	CHECK2=$?',"\n";

	print STCH '	if [ ${CHECK1} -ne 0 -o ${CHECK2} -ne 0  ]',"\n"; # value is 0 if completed correctly.
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to process split for BLASTX VIRUSDB.\"  \n";
	print STCH '		exit 1',"\n";
	print STCH "	fi\n";
	print STCH "	touch $status_log/j15_finished_split_for_BLASTX_VIRUSDB \n";
	print STCH "fi\n";
	close STCH;

	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 16: Run blastx against viral-only database##########
sub submit_job_array_BLASTX_VIRUSDB{
	# this is the job script
	my $job_script = $job_files_dir."/j16_".$sample_name."_BLASTX_VIRUSDB_run_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "#!/bin/bash\n";
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";        
	print STCH "VIRUS_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "BLASTX_VIRUSDB_OUT=",'${VIRUS_DIR}',"/",$sample_name.".MMseqs_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BLASTX_VIRUSDB.out\n";#full path
	print STCH "QUERY=",'${VIRUS_DIR}',"/".$sample_name.".MMseqs_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n\n";
         
	print STCH 'if [ -s $QUERY ]',"\n"; # check if the file is empty
	print STCH "then\n";
	print STCH "	if [ ! -e $status_log/j16_\${SLURM_ARRAY_TASK_ID}_finished_BLASTX_VIRUSDB ]\n";
	print STCH "	then\n";
	print STCH "		$diamond --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen --threads \$SLURM_CPUS_PER_TASK --evalue 1e-2 --query \${QUERY} --out \${BLASTX_VIRUSDB_OUT} --db $db_NR_VIRUS ","\n";
	print STCH '		CHECK1=$?',"\n";
	print STCH '		if [ ${CHECK1} -ne 0 ]',"\n"; 
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to run diamond blastx VIRUSDB.\"  \n";
	print STCH '			exit 1',"\n";
	print STCH '		else', "\n"; 
	print STCH "			touch $status_log/j16_\${SLURM_ARRAY_TASK_ID}_finished_BLASTX_VIRUSDB \n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;
	
 	# here to submit
	$current_job_file1 = "j16_".$sample_name."_BLASTX_VIRUSDB.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage16\n";
	print STCH "#SBATCH --mem=${slurm_mem}\n";
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX_VIRUSDB\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	#print STCH "set -x\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 17: Parse output of blastx against viral-only database##########
sub parse_BLASTX_VIRUSDB {
	$current_job_file1 = "j17_".$sample_name."_Parse_BLASTX_VIRUSDB.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage17\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX_VIRUSDB\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	#print STCH "set -x\n";
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "VIRUS_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."/\n";
	print STCH "BLASTX_VIRUSDB_OUT=",$sample_name.".MMseqs_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BLASTX_VIRUSDB.out\n"; #name only, not full path
	print STCH "BLASTX_VIRUSDB_IN=",$sample_name.".MMseqs_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n"; #full path
	print STCH "PARSED=",'${VIRUS_DIR}',"/".$sample_name.".MMseqs_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BLASTX_VIRUSDB.parsed\n\n";
	print STCH "BLASTX_VIRUSDB_Filtered_fa=",'${VIRUS_DIR}',"/".$sample_name.".MMseqs_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BXVIRUSDB_filtered.fa\n\n";       
	
	print STCH 'if [ -f ${VIRUS_DIR}/${BLASTX_VIRUSDB_OUT} ]',"\n"; # input file exist
	print STCH "then\n";
	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -f $PARSED ]',"\n";
	print STCH "	then\n";
        print STCH "            python ".$bash_bin_path."/parse_blast_Virus_noBP.py -d "."\${VIRUS_DIR} -i \${BLASTX_VIRUSDB_OUT} -f \${BLASTX_VIRUSDB_IN}\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl  \${PARSED}\n"; 
	print STCH '		CHECK=$?',"\n"; # if finished correctly the value should be 0.
	
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # parsing did not complete correctly  
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${BLASTX_VIRUSDB_OUT}.\"  \n";
	print STCH "			exit 1\n"; # exit
	print STCH "		fi\n";
	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # parsing did not complete correctly  
	print STCH "		then\n";
	print STCH "			python ".$bash_bin_path."/parse_blast_Virus_noBP.py -d "."\${VIRUS_DIR} -i \${BLASTX_VIRUSDB_OUT} -f \${BLASTX_VIRUSDB_IN}\n";
	print STCH "			".$run_script_path."check_Blast_parsed_file.pl \${PARSED}\n";
	print STCH '			CHECK=$?',"\n"; # if finished correctly the value should be 0.
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n"; # parsing did not complete correctly  
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${BLASTX_VIRUSDB_OUT}.\"  \n";
	print STCH "				exit 1\n"; # exit
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 18: Split input file for MMseqs against NT##########
sub split_for_MMseqs_NT{
 	$current_job_file1 = "j18_".$sample_name."_MMseqs_NT_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage18\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "MMseqs_DIR=".$MMseqs_dir_NT."\n";
	print STCH "VIRUSDB_HIT_all=".$sample_name.".VIRUSDB_HIT.fa\n";
	print STCH "MMseqs_VIRUSDB_Hit_fa=".$sample_name.".MMseqs_VIRUSDB_HIT.fa\n";
	print STCH "BLASTX_VIRUSDB_Hit_fa=".$sample_name.".BLASTX_VIRUSDB_HIT.fa\n";
	print STCH "BLASTX_VIRUSDB_Phage_fa=".$sample_name.".BLASTX_VIRUSDB_Phage.fa\n";
	print STCH "BLASTX_VIRUSDB_Filtered_fa=".$sample_name.".BLASTX_VIRUSDB_Filtered.fa\n";
	print STCH "BLASTX_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."\n\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".BLASTX_VIRUSDB_finish_check.txt\n\n";


	print STCH "if [ ! -e $status_log/j18_finished_MMseqs_NT_split ]\n"; # the job never finished
	print STCH "then\n";

	# check to make sure the number of reads add together equals the total input reads
	print STCH "	".$run_script_path."check_BLASTX_VIRUSDB_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : BLASTX_VIRUSDB did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	# pool all the BLASTX VIRUSDB hits into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} \n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTX_VIRUSDB_DIR}/*.BXVIRUSDB_hit.fa >> $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} \n\n";

	# pool all the BLASTX VIRUSDB filtered into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTX_VIRUSDB_Filtered_fa} ] \n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTX_VIRUSDB_Filtered_fa} \n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTX_VIRUSDB_DIR}/*.BXVIRUSDB_filtered.fa >> $sample_dir/\${BLASTX_VIRUSDB_Filtered_fa} \n\n";

	# pool all the BLASTX VIRUSDB phage hits into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTX_VIRUSDB_Phage_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTX_VIRUSDB_Phage_fa} \n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTX_VIRUSDB_DIR}/*.BXVIRUSDB_phage.fa >> $sample_dir/\${BLASTX_VIRUSDB_Phage_fa} \n\n";

	# pool BLASTX VIRUSDB hits and MMseqs VIRUSDB hits into a single file
	print STCH "	if [ -e $sample_dir/\${VIRUSDB_HIT_all} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${VIRUSDB_HIT_all} \n";
	print STCH "	fi\n";
	print STCH "	cat  $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} $sample_dir/\${MMseqs_VIRUSDB_Hit_fa} > \${SAMPLE_DIR}/\${VIRUSDB_HIT_all} \n";


	# if the directory already exist, remove it.
	print STCH "	if [ -d \${MMseqs_DIR} ]\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${MMseqs_DIR}\n";
	print STCH "	fi\n";

	# make the directory
	print STCH "	mkdir -p \${MMseqs_DIR}\n\n";

	# split into smaller files
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${VIRUSDB_HIT_all} -o \${MMseqs_DIR}   -f $file_number_of_MMseqs_NT -n $num_seq_per_file_MMseqs_NT \n";
	print STCH "	".$run_script_path."check_split_MMseqs_NT.pl $sample_dir \n";
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : split_for_MMseqs_NT did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	touch $status_log/j18_finished_MMseqs_NT_split \n";
	print STCH "fi\n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 19: Run MMseqs against NT##########
sub submit_job_array_MMseqs_NT {
 	# this is the job script
	my $job_script = $job_files_dir."/j19_".$sample_name."_MMseqs_NT_run_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "#!/bin/bash\n";         
        print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
        print STCH "conda activate vs\n";  
	print STCH "MMseqs_DIR=".$MMseqs_dir_NT."\n";
	print STCH "QUERY=",'${MMseqs_DIR}',"/".$sample_name.".VIRUSDB_HIT.fa_file".'${SLURM_ARRAY_TASK_ID}'.".fasta\n";
	print STCH "MMseqsOUT=",'${MMseqs_DIR}',"/".$sample_name.".VIRUSDB_HIT.fa_file".'${SLURM_ARRAY_TASK_ID}',".mmseqs.out\n";
        print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".mmseqs_j19_time.txt\n";

	print STCH 'if [ -s $QUERY ]',"\n"; 
	print STCH "then\n"; 
	#if mmseqs output file does not exist, do mmseqs 
	print STCH '	if [ ! -f $MMseqsOUT ]',"\n";
	print STCH "	then\n"; 
	print STCH "		date > \${TIMEFILE}\n";
	print STCH "		$mmseqs_easysearch \${QUERY} $db_MMSEQS_NT \${MMseqsOUT} $mmseqs_tmpdir --max-accept 10 --threads \$SLURM_CPUS_PER_TASK --start-sens 5 --sens-steps 3 -s 7 -e 1.0E-4 --split-memory-limit \${SLURM_MEM_PER_NODE} --search-type 3 --format-output \"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen\"","\n";
 	print STCH "		OUT=\$?\n"; 
	print STCH '		if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error : MMseqs \${QUERY}  did not finish correctly.\"  \n";
	print STCH "			exit 1\n"; 
	print STCH "		else\n";
	print STCH "            	date >> \${TIMEFILE}\n";
	print STCH "				#echo \"MMseqs_s5_COMPLETED\" >> \${MMseqsOUT} \n";
 	print STCH "				touch $status_log/j19_MMseqs_NT_run_script_finished \n";
	print STCH "		fi\n";
	#if mmseqs output file exists, check the completeness of output
	print STCH "	else\n";
	print STCH "		if [ ! -f $status_log/j19_MMseqs_NT_run_script_finished ]","\n";
	print STCH "		then\n";
	print STCH "			date > \${TIMEFILE}\n";
	print STCH "			$mmseqs_easysearch \${QUERY} $db_MMSEQS_NT \${MMseqsOUT} $mmseqs_tmpdir --max-accept 10 --threads \$SLURM_CPUS_PER_TASK --start-sens 5 --sens-steps 3 -s 7 -e 1.0E-4 --split-memory-limit \${SLURM_MEM_PER_NODE} --search-type 3 --format-output \"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen\"","\n";
	print STCH "			OUT=\$?\n"; 
	print STCH '			if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error : MMseqs \${QUERY} did not finish correctly.\"  \n";
	print STCH "				exit 1\n";
	print STCH "			else\n";
	print STCH "				date >> \${TIMEFILE}\n";
	print STCH "				touch $status_log/j19_MMseqs_NT_run_script_finished \n";
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

 	# here to submit
	$current_job_file1 = "j19_".$sample_name."_MMseqs_NT.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage19\n";
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --mem=${slurm_mem}\n";
	print STCH "#SBATCH --array=1-$file_number_of_MMseqs_NT\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	#print STCH "set -x\n";

	if ($use_checkpoint) { 
		print STCH "checkpoint bash $job_script\n";
	}
	else {
		print STCH "bash $job_script\n";
	}
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 20: Parse output for MMseqs against NT##########
sub parse_MMseqs_NT{
	$current_job_file1 = "j20_".$sample_name."_Parse_MMseqs_NT.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage20\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "#SBATCH --array=1-$file_number_of_MMseqs_NT\n";
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "MMseqs_DIR=".$MMseqs_dir_NT."\n";
	print STCH "QUERY=".$sample_name.".VIRUSDB_HIT.fa_file".'${SLURM_ARRAY_TASK_ID}'.".fasta\n";
	print STCH "MMseqsOUT=".$sample_name.".VIRUSDB_HIT.fa_file".'${SLURM_ARRAY_TASK_ID}',".mmseqs.out\n"; #name only, not full path
	print STCH "PARSED=\${MMseqs_DIR}/".$sample_name.".VIRUSDB_HIT.fa_file".'${SLURM_ARRAY_TASK_ID}'.".mmseqs.parsed\n\n";
	print STCH "MMseqs_Filtered_fa=".$sample_name.".VIRUSDB_HIT.fa_file".'${SLURM_ARRAY_TASK_ID}',".MMseqs_filtered.fa\n\n";

	print STCH 'if [ -f ${MMseqs_DIR}/${MMseqsOUT} ]',"\n"; 
	print STCH "then\n";

	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -f $PARSED ]',"\n"; # no parsed file
	print STCH "	then\n";
	print STCH '		python '.$bash_bin_path."/parse_blast_NTNR_noBP.py -d \${MMseqs_DIR}/ -f \${QUERY} -i \${MMseqsOUT}\n";
	print STCH "		".$run_script_path."check_MMseqs_parsed_file.pl \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${MMseqsOUT}.\"  \n";
	print STCH "			exit 1\n";
	print STCH "   		fi\n";

	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_MMseqs_parsed_file.pl  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, rerun parser.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH '			python '.$bash_bin_path."/parse_blast_NTNR_noBP.py -d \${MMseqs_DIR}/ -f \${QUERY} -i \${MMseqsOUT}\n";
	print STCH "			".$run_script_path."check_MMseqs_parsed_file.pl  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${MMseqsOUT}.\"  \n";
	print STCH "				exit 1\n";
	print STCH "   			fi\n";
	print STCH "   		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 21: Split input file for diamond blastx against NR##########
sub pool_split_for_BLASTX{
	$current_job_file1 = "j21_".$sample_name."_BLASTX_NR_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage21\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "BX_DIR=".$BLASTX_NR_dir."\n";
	print STCH "MMSeqs_DIR=".$MMseqs_dir_NT."\n";
	print STCH "MMSeqs_Filtered_fa=".$sample_name.".MMseqs_NT_Filtered.fa\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".MMseqs_NT_check.txt\n\n";
	print STCH "NO_MORE_SEQS=",$status_log."/j21_MMseqs_filtered_file_is_empty\n\n";

	print STCH 'if [ -e $NO_MORE_SEQS ]',"\n";
	print STCH "then\n";
	print STCH "    touch $status_log/j21_BX_split_finished  \n";
	print STCH '    mkdir -p ${BX_DIR}',"\n";
	print STCH "    exit\n";
	print STCH "fi\n\n";

	print STCH "if [ ! -e $status_log/j21_BX_split_finished ]\n"; # job never finished
	print STCH "then\n";
	# check to make sure all parser finished 
	print STCH "	".$run_script_path."check_MMseqs_NT_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : MMseqs did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	print STCH "	if [ -d \${BX_DIR} ]\n"; # If the dir already exists, remove.
	print STCH "	then\n";
	print STCH "		rm -rf \${BX_DIR}\n";
	print STCH "	fi\n\n";
	print STCH "	mkdir -p \${BX_DIR}\n";

	print STCH "	if [ -e \${MMSeqs_Filtered_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm \$MMSeqs_Filtered_fa\n";
	print STCH "	fi\n\n";
	print STCH "	cat $MMseqs_NT_dir/*.MMseqs_filtered.fa > $sample_dir/\${MMSeqs_Filtered_fa}\n";

	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${MMSeqs_Filtered_fa}  -o \${BX_DIR} -f $file_number_of_BLASTX -n $num_seq_per_file_BLASTX \n";
	print STCH "	".$run_script_path."check_split_BLASTX_NR.pl ".$sample_dir."\n";
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : change unmasked sequence to masked did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	touch $status_log/j21_BX_split_finished  \n";
	print STCH "fi\n";
	 
	#check if filtered file contains anything; if not, end blasts
	print STCH "if [ ! -s \${MMSeqs_DIR}/*MMseqs_filtered.fa ]\n";
	print STCH "then\n";
	print STCH "    touch $status_log/j21_MMSeqs_filtered_file_is_empty\n";
	print STCH "fi";
	close STCH;
	
    $last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 22: Run diamond blastx against NR##########
sub submit_job_array_BLASTX{
	# this is the job script
	my $job_script = $job_files_dir."/j22_".$sample_name."_BLASTX_NR_run_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "#!/bin/bash\n";
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";        
	print STCH "set -x\n";
	print STCH "BX_DIR=".$BLASTX_NR_dir."\n";
	print STCH "BlastXOUT=",'${BX_DIR}',"/",$sample_name.".MMseqs_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastx.out\n";#full path
	print STCH "QUERY=",'${BX_DIR}',"/".$sample_name.".MMseqs_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n\n";
	print STCH "NO_MORE_SEQS=",$status_log."/j21_MMSeqs_filtered_file_is_empty\n\n";

	print STCH 'if [ -e $NO_MORE_SEQS ]',"\n";
	print STCH "then\n";
	print STCH "    exit\n";
	print STCH "fi\n\n";
	
	print STCH 'if [ -s $QUERY ]',"\n"; 
	print STCH "then\n";
	print STCH "	if [ ! -e $status_log/j22_\${SLURM_ARRAY_TASK_ID}_finished_BLASTX_NR ]\n";
	print STCH "	then\n";
	print STCH "		$diamond --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen --threads \$SLURM_CPUS_PER_TASK --evalue 1e-4 --query \${QUERY} --out \${BlastXOUT} --db $db_NR ","\n";
	print STCH '		CHECK1=$?',"\n";
	print STCH '		if [ ${CHECK1} -ne 0 ]',"\n"; 
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error : BLASTx \${QUERY} did not finish correctly.\"  \n";
	print STCH "			exit 1\n"; 
	print STCH '		else',"\n"; 
	print STCH "			touch $status_log/j22_\${SLURM_ARRAY_TASK_ID}_finished_BLASTX_NR \n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;                                                                                              
  
	# here to submit
	$current_job_file1 = "j22_".$sample_name."_BLASTX_NR.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage22\n";
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --mem=${slurm_mem}\n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	#print STCH "set -x\n";

	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 23: Parse output for blastx against NR##########
sub parse_BLASTX{
	$current_job_file1 = "j23_".$sample_name."_Parse_BLASTX_NR.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage23\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}

	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "BX_DIR=".$BLASTX_NR_dir."\n";
	print STCH "BlastXOUT=",$sample_name.".MMseqs_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastx.out\n";#name only, not full path
	print STCH "BlastXIN=".$sample_name.".MMseqs_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n";
	print STCH "PARSED=",'${BX_DIR}',"/".$sample_name.".MMseqs_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastx.parsed\n";
	print STCH "NO_MORE_SEQS=",$status_log."/j21_MMseqs_filtered_file_is_empty\n\n";
	print STCH "BLASTX_Filtered_fa=",'${BX_DIR}',"/".$sample_name.".MMseqs_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BXfiltered.fa\n\n";

	print STCH 'if [ -e $NO_MORE_SEQS ]',"\n";
	print STCH "then\n";
	print STCH "    exit\n";
	print STCH "fi\n\n";

	print STCH 'if [ -f ${BX_DIR}/${BlastXOUT} ] ',"\n";  
	print STCH "then\n";
	print STCH '	if [ ! -e $PARSED ] ',"\n"; #blastx.parsed file not exist
	print STCH "	then\n";
	print STCH '           python '.$bash_bin_path."/parse_blast_NTNR_noBP.py -d \${BX_DIR} -f \${BlastXIN} -i \${BlastXOUT}\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${BlastXOUT}.\"  \n";
	print STCH "			exit 1\n";
	print STCH "   		fi\n";
	
	print STCH "	else\n"; #blastx.parsed file exists
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl  \${PARSED}\n"; #check if the blastx.parsed file completed
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH '			python '.$bash_bin_path."/parse_blast_NTNR_noBP.py -d \${BX_DIR} -f \${BlastXIN} -i \${BlastXOUT}\n";
	print STCH "			".$run_script_path."check_Blast_parsed_file.pl  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${BlastXOUT}.\"  \n";
	print STCH "				exit 1\n";
	print STCH "   			fi\n";
	print STCH "   		fi\n";
	print STCH "	fi\n";
	print STCH "else\n";
	print STCH '    if [ -s ${BX_DIR}/${BlastXIN} ] ',"\n"; #blastx fasta input file is not empty
	print STCH "    then\n";
        print STCH "            echo what\n";
	#print STCH "            cp \${BX_DIR}/\${BlastXIN} \${BLASTX_Filtered_fa}\n";
	print STCH "    fi\n";
	print STCH "fi";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 24: map reads to viral reference genome##########
sub map_to_viral_ref{
	# this is the job script
	my $job_script = $job_files_dir."/j24_".$sample_name."_map_to_viral_ref_script.sh";
	open(JSTCH, ">$job_script") or die $!;
	chmod 0755, $job_script;
	print JSTCH "#!/bin/bash\n";
	print JSTCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.bbduk.fastq.gz\n";
	print JSTCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print JSTCH "conda activate vs\n";

	# Check if mapping finished successfully. If true, resubmitting will skip this step and keep results
	print JSTCH "if [ ! -e $status_log/j24_map_to_viral_ref_finished ]\n";
	print JSTCH "then\n";
	print JSTCH "   date > \${TIMEFILE}\n";
	print JSTCH "   bbsplit.sh in=\${IN} ref=$viral_reference_genome covstats=$sample_dir"."/".$sample_name.".bbsplit_viral_covstats.txt tossbrokenreads fast=t nodisk=t ow pigz=t -Xmx${slurm_mem}\n";
	print JSTCH "   OUT=\$?\n";
	print JSTCH '   if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print JSTCH "   then\n";
	print JSTCH "           echo \"Fatal Error trying to run bbsplit.\"  \n";
	print JSTCH "           exit 1\n";
	print JSTCH "   else\n";
	print JSTCH "           touch $status_log/j24_map_to_viral_ref_finished\n";
	print JSTCH "   fi\n";
	print JSTCH "fi\n\n";
	close JSTCH;

	# here to submit
	$current_job_file1 = "j24_".$sample_name."_map_to_viral_ref.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage24\n";
	print STCH "#SBATCH --mem=${slurm_mem}\n"; 
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
			print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	if ($use_checkpoint) {
			print STCH "checkpoint bash  $job_script \n";
	}
	else {
			print STCH " bash $job_script \n";
	}
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 24B: map reads to viral reference genome##########
sub map_to_viral_ref_b{
	# this is the job script
	my $job_script = $job_files_dir."/j24b_".$sample_name."_map_to_viral_ref_script.sh";
	open(JSTCH, ">$job_script") or die $!;
	chmod 0755, $job_script;
	print JSTCH "#!/bin/bash\n";
	print JSTCH "IN2=$sample_dir"."/".$sample_name.".RemoveAdapter.LR.fastp.fastq.gz\n";
	print JSTCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.LR.fastp.fasta\n";
	print JSTCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print JSTCH "conda activate vs\n";

	# Check if mapping finished successfully. If true, resubmitting will skip this step and keep results
	print JSTCH "if [ ! -e $status_log/j24b_map_to_viral_ref_finished ]\n";
	print JSTCH "then\n";
	print JSTCH "   date > \${TIMEFILE}\n";
	print JSTCH "   seqkit fq2fa \${IN2} > \${IN} \n"; 
	print JSTCH "   bbsplit.sh in=\${IN} ref=$viral_reference_genome covstats=$sample_dir"."/".$sample_name.".bbsplit_viral_covstats.txt tossbrokenreads fast=t nodisk=t ow pigz=t -Xmx${slurm_mem} maxlen=50000\n";
	print JSTCH "   OUT=\$?\n";
	print JSTCH '   if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print JSTCH "   then\n";
	print JSTCH "           echo \"Fatal Error trying to run bbsplit.\"  \n";
	print JSTCH "           exit 1\n";
	print JSTCH "   else\n";
	print JSTCH "           touch $status_log/j24b_map_to_viral_ref_finished\n";
	print JSTCH "   fi\n";
	print JSTCH "fi\n\n";
	close JSTCH;

	# here to submit
	$current_job_file1 = "j24b_".$sample_name."_map_to_viral_ref.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage24b\n";
	print STCH "#SBATCH --mem=${slurm_mem}\n"; 
	print STCH "#SBATCH --cpus-per-task=$num_CPU\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
			print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	if ($use_checkpoint) {
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash $job_script \n";
	}
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}


##########STAGE 25: Generate Assignment Report##########
sub generate_assignment_report{
	$current_job_file1 = "j25_".$sample_name."_generate_report.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage25\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "REPORT=".$sample_dir."/".$sample_name.".AssignmentReport\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".AssignmentReport_time.txt\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".BLASTX_NR_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j25_AssignmentReport_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";

	print STCH "	".$run_script_path."Generate_assignment_report.pl ".$sample_dir." ".$VStrack." \${INPUT}\n";
	print STCH '	grep "# Finished Assignment Report" ${REPORT}',"\n";
	print STCH '	CHECK=$?',"\n";
	print STCH '	if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error in generating report for \${SAMPLE_DIR}.\"  \n";
	print STCH "		exit 1\n";
	print STCH "   	fi\n";
	print STCH "fi\n";
	print STCH "touch $status_log/j25_AssignmentReport_finished  \n";
	print STCH "date >> \${TIMEFILE}\n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 26: Generate Assignment Summary##########
sub generate_assignment_summary {
	$current_job_file1 = "j26_".$sample_name."_generate_summary.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage26\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  

	print STCH "OUTPUT=".$sample_dir."/".$sample_name.".AssignmentSummary\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".AssignmentSummary_time.txt\n";

	print STCH "if [ ! -e $status_log/j26_AssignmentSummary_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";
	print STCH "	".$run_script_path."Generate_assignment_summary.pl ".$sample_dir." \${BAD_SEQ}\n";
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : Generate_assignment_summary did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "fi\n";
	print STCH "date >> \${TIMEFILE}\n";
	print STCH "touch $status_log/j26_AssignmentSummary_finished  \n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 27: Generate Phage Report##########
sub generate_phage_report{
	$current_job_file1 = "j27_".$sample_name."_generate_phage_report.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage27\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";

	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "REPORT=".$sample_dir."/".$sample_name.".PhageAssignmentReport\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".PhageAssignmentReport_time.txt\n";

	print STCH "if [ ! -e $status_log/j27_PhageAssignmentReport_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";
	print STCH "	".$run_script_path."Phage_Generate_assignment_report.pl ".$sample_dir."\n";
	print STCH '	grep "# Finished Assignment Report" ${REPORT}',"\n";
	print STCH '	CHECK=$?',"\n";
	print STCH '	if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error in generating phage report for \${SAMPLE_DIR}.\"  \n";
	print STCH "		exit 1\n";
	print STCH "   	fi\n";
	print STCH "fi\n";
	print STCH "touch $status_log/j27_PhageAssignmentReport_finished \n";
	print STCH "date >> \${TIMEFILE}\n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 28: Generate Phage Summary##########
sub generate_phage_summary {
	$current_job_file1 = "j28_".$sample_name."_generate_phage_summary.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage28\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "OUTPUT=".$sample_dir."/".$sample_name.".PhageAssignmentSummary\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".PhageAssignmentSummary_time.txt\n";
	print STCH "if [ ! -e $status_log/j28_PhageAssignmentSummary_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";
	print STCH '	CHECK=1',"\n";
	print STCH '	while [ ${CHECK} -ne 0 ] ',"\n"; 
	print STCH "	do\n";
	print STCH "		".$run_script_path."Phage_Generate_assignment_summary.pl ".$sample_dir."\n";
	print STCH '		grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
	print STCH '		CHECK=$?',"\n";
	print STCH "	done\n";
	print STCH "fi\n";
	print STCH "date >> \${TIMEFILE}\n";
	print STCH "	touch $status_log/j28_PhageAssignmentSummary_finished  \n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##########STAGE 29: Generate Supplemental Outputs##########
sub generate_supplemental_read_counts {
	$current_job_file1 = "j29_".$sample_name."_VS_Read_Counter.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --partition=$partition\n";
	print STCH "#SBATCH --account=$account\n";
	print STCH "#SBATCH --job-name=${sample_path}:stage29\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number or $step_number == -1) {
		print STCH "#SBATCH --depend=afterok:$last_jobid1\n";
	}
	print STCH "source ${conda_prefix}/etc/profile.d/conda.sh\n";
	print STCH "conda activate vs\n";  
	print STCH "if [ ! -e $status_log/j29_VS_Read_Counter_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	".$bash_bin_path."/VS_Read_Counter/run_VS_Read_Counter.sh ".$sample_dir." ".$bash_bin_path."\n";
	print STCH '	CHECK=$?',"\n";
	print STCH '	if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "	then\n";
	print STCH "		echo Fatal Error trying to run VS_Read_Counter.sh  \n";
	print STCH "		exit 1\n";
	print STCH "	else\n";
	print STCH "		touch $status_log/j29_VS_Read_Counter_finished  \n";
	print STCH "   	fi\n";
	print STCH "fi\n";
	close STCH;
	$last_jobid1 = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}