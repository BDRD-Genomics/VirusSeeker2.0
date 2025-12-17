###############################################################################
# Author: BDRD <usn.detrick.nmrc.mbx.genomics-reach-back@health.mil>
#
# License:
# VirusSeeker 2.0 is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3 of the License or any later version. For 
# details, please refer to https://www.gnu.org/licenses/
###############################################################################



#!/usr/bin/perl
use strict;
use Switch;
use Bio::SearchIO;

my $usage = '
This script will read corresponding files in the given director and 
generate a report. 

perl script <sample dir> <VS track>
<sample dir> = full path to the directory holding files for the given 
               library without the last "/"
               e.g. ~/tools/Illumina_VirusDiscoveryPipeline/data/MiSeq_run_1/I10_12310_Project404_CSF_Glaser_Encephalitis_V11T00919_TruSeq-AD007_CAGATC
<VS track> = Discovery (D) or Virome (V)
';

die $usage unless scalar @ARGV == 2;
my ( $dir, $VStrack ) = @ARGV;

my @temp = split("\/", $dir);
my $sample_name = pop @temp;
# print "lib is $sample_name\n";

####################################
# read in original sequences
my $fasta_file = $dir."/".$sample_name.".allSeqs.mmseqs_rep_seq.fasta"; #
my $seq_dir = $dir."/".$sample_name."_all_seqs";
my %seq = (); # read_ID => sequence
my %seq_desc = (); # read_ID => description
&read_FASTA_data($fasta_file, \%seq, \%seq_desc);

####################################
# output files 
my $out1 = $dir."/".$sample_name.".AssignmentReport";
open (OUT1, ">$out1") or die "can not open file $out1!\n";
my $out2 = $dir."/".$sample_name.".ViralReads_all.fa";
open (OUT2, ">$out2") or die "can not open file $out2!\n";
my $out3 = $dir."/".$sample_name.".AmbiguousReadsAssignmentReport";
open (OUT3, ">$out3") or die "can not open file $out3!\n";
my $out4 = $dir."/".$sample_name.".AmbiguousReads_all.fa";
open (OUT4, ">$out4") or die "can not open file $out4!\n";

##############################################################################
# variable to keep all the blast output files
my @blast_files_MMseqs_NT = (); # all MMseqs output files
#my @blast_files_blastn = (); # all blastn.out files
my @blast_files_blastx = (); # all blastx.out files

# all the viral read ID
my %viral_reads_MMseqs_NT = ();
my %viral_reads_blastx = ();
my %viral_reads_blastx_VIRUSDB = ();

# viral_lineage => number of reads assigned to this lineage in the sample
my %viral_lineage_to_num_reads = ();
my %viral_reads_blast_info =();    # readID => information about this read
my %viral_reads_in_lineage_MMseqs_NT = ();    # lineage => [read ID]
my %viral_reads_in_lineage_blastx = ();    # lineage => [read ID]

# all the ambiguous viral read ID
my %ambiguous_reads_MMseqs_NT = ();
my %ambiguous_reads_blastx = ();

# ambiguous reads viral lineage => number of reads assigned to this lineage in the sample
my %ambiguous_reads_lineage_to_num_reads = ();
my %ambiguous_reads_blast_info =();    # readID => information about this read
my %ambiguous_reads_in_lineage_MMseqs_NT = ();    # lineage => [read ID]
my %ambiguous_reads_in_lineage_blastx = ();    # lineage => [read ID]

my @unassigned_reads = ();
my @phage_reads = ();
##################################
# generate unassigned reads file
my $BLASTX_NR_dir = $dir."/".$sample_name."_BLASTX_NR";
my $unassigned_reads_file = $dir."/".$sample_name.".unassigned_reads.fa";
if (-f $unassigned_reads_file) {
	unlink $unassigned_reads_file;
}
my $com = "cat $BLASTX_NR_dir/*BXfiltered.fa >> $unassigned_reads_file";
system ($com);


# category => num of sequence assigned to this category template
my %Assignment_temp = ( 
	"Bacteria" => 0,
	"Fungi" => 0,
	"Homo" => 0,
	"Mus" => 0,
	"Phage" => 0,
	"Viruses" => 0,
	"other" => 0,
	"unassigned" => 0,
	"Ambiguous" => 0,
);

# category => num of sequence assigned to this category by MMseqs 
my %Assignment_MMseqs_NT = ();
foreach my $key (keys %Assignment_temp) {
	$Assignment_MMseqs_NT{$key} = 0;
}
# category => num of sequence assigned to this category by blastx 
my %Assignment_blastx = ();
foreach my $key (keys %Assignment_temp) {
	$Assignment_blastx{$key} = 0;
}

# category => num of sequence assigned to this category by blastx VIRUSDB 
my %Assignment_blastx_VIRUSDB = ();
foreach my $key (keys %Assignment_temp) {
	$Assignment_blastx_VIRUSDB{$key} = 0;
}

######################################################################################
# get statistics of the sample
##########Total raw number of reads##########
my $SE1 = $dir."/".$sample_name."_SE1.RemoveAdapter.fastq.gz";
my $line_number_SE1 = `zcat $SE1 | wc -l`;
my $SE1_total_NumReads_raw = $line_number_SE1/4;
my $SE2 = $dir."/".$sample_name."_SE2.RemoveAdapter.fastq.gz";
my $line_number_SE2 = `zcat $SE1 | wc -l`;
my $SE2_total_NumReads_raw = $line_number_SE2/4;
my $SE_total_NumReads_raw = $SE1_total_NumReads_raw + $SE2_total_NumReads_raw;

##########Number of reads after quality control##########
my $file = $dir."/".$sample_name.".RemoveAdapter.bbduk.fastq.gz";
my $line_number = `zcat $file | wc -l`;
my $NumReads_QC = $line_number/4;

##########Number of reads not mapped to reference genome##########
$SE1 = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.fastq.gz";
my $line_number_unmapped = `zcat $SE1 | wc -l`;
my $NumReads_unmapped = $line_number_unmapped/4;

##########Number of reads mapped to reference genome########## 
$SE1 = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.mapped.fastq.gz";
my $line_number_mapped = `zcat $SE1 | wc -l`;
my $NumReads_mapped = $line_number_mapped/4;

##########Number of contigs obtained##########
my $file = $seq_dir."/contigs.fasta";
my $num_contigs = `grep ">" $file |wc -l`;

##########Number of reads mapped to assembly##########
my $pe_map = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.mapped_pe.fastq.gz";
my $line_number_mapped_pe = `zcat $pe_map | wc -l`;
my $se_map = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.mapped_se.fastq";
my $line_number_mapped_se = `zcat $se_map | wc -l`;
my $assembly_mapped = $line_number_mapped_pe/4 + $line_number_mapped_se/4;

##########Number of reads not mapped to assembly########## 
my $pe_unmap = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.fastq.gz";
my $line_number_unmapped_pe = `zcat $pe_unmap | wc -l`;
my $se_unmap = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.fastq.gz";
my $line_number_unmapped_se = `zcat $se_unmap | wc -l`;
my $assembly_unmapped = $line_number_unmapped_pe/4 + $line_number_unmapped_se/4;


##########Number of reads stitched##########
my $unstitched_SE1 = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.unassembled.forward.fastq";
my $unstitched_SE2 = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.unassembled.reverse.fastq";
my $joined = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.assembled.fastq";
if ($VStrack eq "V") {
    my $unstitched_SE1 = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped_pe.unassembled.forward.fastq";
    my $unstitched_SE2 = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped_pe.unassembled.reverse.fastq";
    my $joined = $dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped_pe.assembled.fastq";
}
my $line_number_SE1 = `wc -l $unstitched_SE1`;
my $un1_NumReads = $line_number_SE1/4;
my $line_number_SE2 = `wc -l $unstitched_SE2`;
my $un2_NumReads = $line_number_SE2/4;
my $line_number_join = `wc -l $joined`;
my $joined_NumReads = $line_number_join/4;


##########Number of singletons##########
my $file = $seq_dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_pe.se.fasta\n";
my $singleton1 = `grep -c ">" $file`;
my $file = $seq_dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.assembly.unmapped_se.fasta\n";
my $singleton2 = `grep -c ">" $file`;
my $total_singleton = $un1_NumReads + $un2_NumReads + $singleton1 + $singleton2;
if ($VStrack eq "V"){
    my $file = $seq_dir."/".$sample_name.".RemoveAdapter.bbduk.RefGenome.unmapped.se.fasta";
    my $singleton1 = `grep -c ">" $file`;
    my $total_singleton = $un1_NumReads + $un2_NumReads + $singleton1;
}

##########Number of sequences before/after clustering##########
$file = $dir."/".$sample_name.".allSeqs.fasta";
my $preClust = `grep -c ">" $file`;
$SE1 = $dir."/".$sample_name.".allSeqs.mmseqs_rep_seq.fasta";
my $postClust = `grep -c ">" $SE1`;

##########Number of good reads after RepeatMasker##########
$file = $dir."/".$sample_name.".goodSeq.fa";
my $RepeatMasker_goodSeq = `grep ">" $file |wc -l`;

##########Number of filtered and low complexity reads after RepeatMasker##########
$file = $dir."/".$sample_name.".badSeq.fa";
my $RepeatMasker_input_total = 0;
my $RepeatMasker_good = 0;
my $RepeatMasker_filtered = 0;
my $RepeatMasker_LowComplexity = 0;

my $temp = `grep "total seq = " $file`;
if ($temp =~ /total seq = (\d+)/) {
	$RepeatMasker_input_total = $1;
}

my $temp = `grep "good seq = " $file`;
if ($temp =~ /good seq = (\d+)/) {
	$RepeatMasker_good = $1;
}

my $temp = `grep "bad seq = " $file`;
if ($temp =~ /bad seq = (\d+)/) {
	$RepeatMasker_filtered = $1;
}
my $temp = `grep "Repeat and Low complexicity seq = " $file`;
if ($temp =~ /Repeat and Low complexicity seq = (\d+)/) {
	$RepeatMasker_LowComplexity = $1;
}

##########Number of viral hits from MMseqs VIRUSDB##########
$file = $dir."/".$sample_name.".MMseqs_VIRUSDB_HIT.fa";
my $MMseqs_VIRUSDB_Hit_number = `grep ">" $file |wc -l`;

##########Number of viral hits from BLASTX VIRUSDB##########
$file = $dir."/".$sample_name.".BLASTX_VIRUSDB_HIT.fa";
my $BLASTX_VIRUSDB_Hit_number = `grep ">" $file |wc -l`;

##########Number of viral hits from MMseqs VIRUSDB and BLASTX VIRUSDB step##########
$file = $dir."/".$sample_name.".VIRUSDB_HIT.fa";
my $VIRUSDB_Hit_number = `grep ">" $file | wc -l`;

##########Number of sequences after MMseqs NT##########
$file = $dir."/".$sample_name.".VIRUSDB_HIT_MMseqs_NT_Filtered.fa";
my $MMseqs_Filtered = `grep ">" $file |wc -l`;

##########Mumber of sequences after BLASTN NT filter and classified at this step##########
#$file = $dir."/".$sample_name.".BLASTN_NT_Filtered.fa";
#my $num_BLASTN_Filtered = `grep ">" $file |wc -l`;
#my $num_BLASTN_Classified = $MMseqs_Filtered - $num_BLASTN_Filtered;

##########Number of sequences after BLASTX NR filter and classified at this step##########
$file = $dir."/".$sample_name.".unassigned_reads.fa";
my $num_unassigned = `grep ">" $file |wc -l`;
my $num_BLASTX_NR_Classified = $MMseqs_Filtered - $num_unassigned;

####################################################################################
# to obtain viral read information
##########MMseqs NT##########
my $MMseqs_dir = $dir."/".$sample_name."_MMseqs_NT";
# enter directory where blast results resides
opendir (DIR, $MMseqs_dir) or die "can not open dir $MMseqs_dir!\n";
foreach my $blast_file (readdir DIR) {
	if ($blast_file =~ /mmseqs\.parsed$/) {
		# print "blastn parsed file $blast_file\n";
		my $parsed = $MMseqs_dir."/".$blast_file;
		my $blast_out = $blast_file;
		$blast_out =~ s/\.mmseqs\.parsed/\.mmseqs\.out/;
		$blast_out = $MMseqs_dir."/".$blast_out;
		push @blast_files_MMseqs_NT, $blast_out;

		&collect_information($parsed, \%Assignment_MMseqs_NT, \%viral_reads_MMseqs_NT, \%viral_reads_in_lineage_MMseqs_NT, \%viral_lineage_to_num_reads, \%ambiguous_reads_MMseqs_NT, \%ambiguous_reads_in_lineage_MMseqs_NT, \%ambiguous_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%ambiguous_reads_blast_info, \@phage_reads );
	}
}
closedir DIR;


##########BLASTX NR##########
my $BLASTX_dir = $dir."/".$sample_name."_BLASTX_NR";
# enter directory where blastx results resides
opendir (BXDIR, $BLASTX_dir) or die "can not open dir $BLASTX_dir!\n";
foreach my $blast_file (readdir BXDIR) {
	if ($blast_file =~ /blastx\.parsed$/) {
		# print "blast parsed file $blast_file\n";
		my $blast_out = $blast_file;
		$blast_out =~ s/\.blastx\.parsed/\.blastx\.out/;
		$blast_out = $BLASTX_dir."/".$blast_out;
		push @blast_files_blastx, $blast_out;
		my $parsed = $BLASTX_dir."/".$blast_file;
		&collect_information($parsed, \%Assignment_blastx, \%viral_reads_blastx, \%viral_reads_in_lineage_blastx, \%viral_lineage_to_num_reads, \%ambiguous_reads_blastx, \%ambiguous_reads_in_lineage_blastx, \%ambiguous_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%ambiguous_reads_blast_info );
	}
}
closedir BXDIR;

=head1
# BLASTX VIRUSDB step
my $BLASTX_VIRUSDB_dir = $dir."/".$sample_name."_BLASTX_VIRUSDB";
# enter directory where blastx results resides
opendir (BXDIR, $BLASTX_VIRUSDB_dir) or die "can not open dir $BLASTX_VIRUSDB_dir!\n";
foreach my $blast_file (readdir BXDIR) {
	if ($blast_file =~ /BLASTX_VIRUSDB\.parsed$/) {
		# print "blast parsed file $blast_file\n";
		my $blast_out = $blast_file;
		$blast_out =~ s/\.BLASTX_VIRUSDB\.parsed/\.BLASTX_VIRUSDB\.out/;
		$blast_out = $BLASTX_VIRUSDB_dir."/".$blast_out;
		push @blast_files_blastx_VIRUSDB, $blast_out;
		my $parsed = $BLASTX_VIRUSDB_dir."/".$blast_file;
		&collect_information($parsed, \%Assignment_blastx_VIRUSDB, \%viral_reads_blastx_VIRUSDB, \%viral_reads_in_lineage_blastx_VIRUSDB, \%viral_lineage_to_num_reads, \%ambiguous_reads_blastx_VIRUSDB, \%ambiguous_reads_in_lineage_blastx_VIRUSDB, \%ambiguous_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%ambiguous_reads_blast_info );
	}
}
closedir BXDIR;
=cut

#####################################################################################
# print out report for this library
print OUT1 $dir, "\n";

printf OUT1 "#####Raw reads, host removal, and QC stats#####\n";
printf OUT1 "%s\t%s\t%s\t%s\n", "#totalReads", "#totalReads_postQC", "#host_mapped", "#host_unmapped";
printf OUT1 "%d\t%d\t%d\t%d\n\n", $SE_total_NumReads_raw, $NumReads_QC, $NumReads_mapped, $NumReads_unmapped;

printf OUT1 "#####Assembly and stitching stats#####\n";
printf OUT1 "%s\t%s\t%s\t%s\t%s\n", "#contigsAssembled", "#assembly_mapped", "#assembly_unmapped", "#stitched_reads", "#singletons";
printf OUT1 "%d\t%d\t%d\t%d\t%d\n\n", $num_contigs, $assembly_mapped, $assembly_unmapped, $joined_NumReads, $total_singleton;

printf OUT1 "#####Clustering and RepeatMasker stats#####\n";
printf OUT1 "%s\t%s\t%s\t%s\n", "#preCluster_Count", "#postCluster_Count", "#RMgood", "#RMfiltered";
printf OUT1 "%d\t%d\t%d\t%d\n\n", $preClust, $postClust, $RepeatMasker_good, ($RepeatMasker_filtered + $RepeatMasker_LowComplexity);

printf OUT1 "#####Viral DB hits#####\n";
printf OUT1 "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#MMseqs_VIRUSDB_Hit", "#BX_VIRUSDB_Hit", "Total_VIRUSDB_Hit","#BNinput", "%VIRUSDB_Hit", "BXinput", "%VIRUSDB_Hit", "Unassigned", "%VIRUSDB_Hit"; 
printf OUT1 "%d\t%7d\t%7d\t%7d\t%7.2f\t%7d\t%7d\t%7.2f\t%7.2f\n\n", $MMseqs_VIRUSDB_Hit_number,  $BLASTX_VIRUSDB_Hit_number, $VIRUSDB_Hit_number, $MMseqs_Filtered, $MMseqs_Filtered*100/$VIRUSDB_Hit_number, $num_unassigned, $num_unassigned*100/$VIRUSDB_Hit_number;

# assignment in different categories
printf OUT1 "%12s\t%7s\t%9s\t%7s\t%7s\t%7s\n", "category", "total", "MMseqs", "BX_NR" ;
foreach my $key (sort {$a cmp $b } keys %Assignment_blastx) {
	printf OUT1 "%12s\t%7d\t%9d\t%7d\t%7d\n", $key, $Assignment_blastx{$key}+$Assignment_MMseqs_NT{$key},  $Assignment_MMseqs_NT{$key},  $Assignment_blastx{$key};
}

##############################################################################
# viral lineage, number of viral reads in this lineage, BLAST alignnment statistics.

my $c1 = "############################################################\n\n";
print OUT1 $c1;

foreach my $lineage (sort {$viral_lineage_to_num_reads{$a} <=> $viral_lineage_to_num_reads{$b}} keys %viral_lineage_to_num_reads) {
	print OUT1 $lineage, "\ttotal number of reads: ", $viral_lineage_to_num_reads{$lineage}, "\n\n";
	print OUT1 "QueryName\tQuerylength\tHitName\tHitLen\tHitDesc \tAlnLen\t%ID\tHitStart\tHitEnd\te\n";

	if (defined $viral_reads_in_lineage_MMseqs_NT{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_MMseqs_NT{$lineage}}) {
			print OUT1 "reads from MMseqs NT:\n";
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_MMseqs_NT{$lineage}}) {
				print OUT1 $viral_reads_blast_info{$read}, "\n";
			}
		}
		print OUT1 "\n";
	}

	if (defined $viral_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastx{$lineage}}) {
			print OUT1 "reads from blastx:\n";
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastx{$lineage}}) {
				print OUT1 $viral_reads_blast_info{$read}, "\n";
			}
		}
		print OUT1 "\n";
	}

=head1
	if (defined $viral_reads_in_lineage_blastx_VIRUSDB{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
			print OUT1 "reads from blastx VIRUSDB:\n";
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
				print OUT1 $viral_reads_blast_info{$read}, "\n";
			}
		}
		print OUT1 "\n";
	}
=cut

	# print out viral sequences in this taxonomy lineage
	if (defined $viral_reads_in_lineage_MMseqs_NT{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_MMseqs_NT{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_MMseqs_NT{$lineage}}) {
				print OUT1 ">$read $seq_desc{$read}\n";
				print OUT1 $seq{$read}, "\n";
			}
		}
	}

	if (defined $viral_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastx{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastx{$lineage}}) {
				print OUT1 ">$read $seq_desc{$read}\n";
				print OUT1 $seq{$read}, "\n";
			}
		}
	}
=head1
	if (defined $viral_reads_in_lineage_blastx_VIRUSDB{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
				print OUT1 ">$read $seq_desc{$read}\n";
				print OUT1 $seq{$read}, "\n";
			}
		}
	}
=cut

	print OUT1 $c1;
}

# generate ambiguous reads report
print OUT3 $c1;
foreach my $lineage (sort {$ambiguous_reads_lineage_to_num_reads{$a} <=> $ambiguous_reads_lineage_to_num_reads{$b}} keys %ambiguous_reads_lineage_to_num_reads) {
	print OUT3 $lineage, "\ttotal number of reads: ", $ambiguous_reads_lineage_to_num_reads{$lineage}, "\n\n";
	print OUT3 "QueryName\tQuerylength\tHitName\tHitLen\tHitDesc\tAlnLen\t%ID\tHitStart\tHitEnd\te\n";

	if (defined $ambiguous_reads_in_lineage_MMseqs_NT{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_MMseqs_NT{$lineage}}) {
			print OUT3 "reads from MMseqs NT:\n";
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_MMseqs_NT{$lineage}}) {
				print OUT3 $ambiguous_reads_blast_info{$read}, "\n";
			}
		}
	}


	if (defined $ambiguous_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
			print OUT3 "reads from blastx:\n";
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
				print OUT3 $ambiguous_reads_blast_info{$read}, "\n";
			}
		}
		print OUT3 "\n";
	}

=head1
	if (defined $ambiguous_reads_in_lineage_blastx_VIRUSDB{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
			print OUT3 "reads from blastx VIRUSDB:\n";
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
				print OUT3 $ambiguous_reads_blast_info{$read}, "\n";
			}
		}
	}
=cut
	##########################################################################
	# print out sequences
	if (defined $ambiguous_reads_in_lineage_MMseqs_NT{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_MMseqs_NT{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_MMseqs_NT{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

	if (defined $ambiguous_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

=head1
	if (defined $ambiguous_reads_in_lineage_blastx_VIRUSDB{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}
=cut
	print OUT3 $c1;
}


###############################################################################
# get all the viral reads and put the sequence into output file:
foreach my $lineage (keys %viral_lineage_to_num_reads) {
	foreach my $read (@{$viral_reads_in_lineage_MMseqs_NT{$lineage}}) {
		print OUT2 ">$read $seq_desc{$read}\n";
		print OUT2 $seq{$read}, "\n";
	}
	# foreach my $read (@{$viral_reads_in_lineage_blastn{$lineage}}) {
	# 	print OUT2 ">$read $seq_desc{$read}\n";
	# 	print OUT2 $seq{$read}, "\n";
	# }
	foreach my $read (@{$viral_reads_in_lineage_blastx{$lineage}}) {
		print OUT2 ">$read $seq_desc{$read}\n";
		print OUT2 $seq{$read}, "\n";
	}
#	foreach my $read (@{$viral_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
#		print OUT2 ">$read $seq_desc{$read}\n";
#		print OUT2 $seq{$read}, "\n";
#	}
}	
	
###############################################################################
# get all the ambiguous reads and put the sequence into output file:
foreach my $lineage (keys %ambiguous_reads_lineage_to_num_reads) {
	foreach my $read (@{$ambiguous_reads_in_lineage_MMseqs_NT{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}

	foreach my $read (@{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}
}

print OUT1 "# Finished Assignment Report\n";

exit;

#####################################################################################
# collect information from given BLAST parsed output file
sub collect_information {
	################################## modified 04/2014
	my ($blast_parsed_file, $category_hash_ref, $viral_reads_hash_ref, $viral_reads_in_lineage_hash_ref, $viral_lineage_to_num_reads_hash_ref,  $ambiguous_reads_hash_ref, $ambiguous_reads_in_lineage_hash_ref, $ambiguous_reads_lineage_to_num_reads_hash_ref, $unassigned_reads_arr_ref, $viral_reads_blast_info_ref, $ambiguous_reads_blast_info_ref ) = @_;

	##################################
	open (IN, $blast_parsed_file) or die "can not open file $blast_parsed_file!\n";
	while (<IN>) {
		if ($_ =~ /^#/) { # skip comment line
			next;
		}
		chomp;
		my @info = split("\t", $_);
		my $read_ID = shift @info;
		my $length = shift @info;
		my $category = shift @info;
		my $lineage = shift @info;
		my $desc = join ("\t", @info);

#		my ($read_ID, $length, $category, $lineage, $hit_name, $hit_length, $hit_desc, $hsp_len, $e_value) = split("\t", $_);
#		print "readID = \\$read_ID\\, length = \\$length\\, category = \\$category\\, lineage = \\$lineage\\, desc = \\$desc\\\n";
		switch ($category ) {
			case "Bacteria" { $category_hash_ref->{"Bacteria"}++	}
			case "Fungi" { $category_hash_ref->{"Fungi"}++ }
			case "Homo" { $category_hash_ref->{"Homo"}++ }
			case "Mus" { $category_hash_ref->{"Mus"}++ }
			case "Phage" {$category_hash_ref->{"Phage"}++ }
			case "Viruses" { $category_hash_ref->{"Viruses"}++ }
			case "other" {$category_hash_ref->{"other"}++ }
			case "unassigned" {$category_hash_ref->{"unassigned"}++}
			case "Ambiguous" {$category_hash_ref->{"Ambiguous"}++ } 
		}

		$desc = $read_ID."\t".$length."\t".$desc; 
		if ($category eq "Viruses") {
			$viral_reads_hash_ref->{$read_ID} = 1;
			$viral_reads_blast_info_ref->{$read_ID} = $desc;
			if (!(defined $viral_reads_in_lineage_hash_ref->{$lineage})) {
				$viral_reads_in_lineage_hash_ref->{$lineage} = [$read_ID];
			}
			else {
				push @{$viral_reads_in_lineage_hash_ref->{$lineage}}, $read_ID;
			}
		
			if (defined $viral_lineage_to_num_reads_hash_ref->{$lineage}) {
				$viral_lineage_to_num_reads_hash_ref->{$lineage}++;
			}
			else {
				$viral_lineage_to_num_reads_hash_ref->{$lineage} = 1;
			}
		################################## modified 2014.04
		}elsif ($category eq "Ambiguous"){
			$ambiguous_reads_hash_ref->{$read_ID} = 1;
			$ambiguous_reads_blast_info_ref->{$read_ID} = $desc;
			if (!(defined $ambiguous_reads_in_lineage_hash_ref->{$lineage})) {
				$ambiguous_reads_in_lineage_hash_ref->{$lineage} = [$read_ID];
			}
			else {
				push @{$ambiguous_reads_in_lineage_hash_ref->{$lineage}}, $read_ID;
			}
		
			if (defined $ambiguous_reads_lineage_to_num_reads_hash_ref->{$lineage}) {
				$ambiguous_reads_lineage_to_num_reads_hash_ref->{$lineage}++;
			}
			else {
				$ambiguous_reads_lineage_to_num_reads_hash_ref->{$lineage} = 1;
			}

		##################################
		}elsif ($category eq "unassigned") {
			push @{$unassigned_reads_arr_ref}, $read_ID;
		}
	}
	close IN;
}

#####################################################################
sub read_FASTA_data () {
    my ($fastaFile, $seq_hash_ref, $seq_desc_hash_ref) = @_;

    #keep old read seperator and set new read seperator to ">"
    my $oldseperator = $/;
    $/ = ">";
	 
    open (FastaFile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";
    while (my $line = <FastaFile>){
		# Discard blank lines
        if ($line =~ /^\s*$/) {
			next;
		}
		# discard the first line which only has ">", keep the rest
		elsif ($line ne ">") {
		    chomp $line;
		    my @rows = ();
		    @rows = split (/\n/m, $line);	
		    my $temp = shift @rows;
			my @temp_arr = split(/\s/, $temp);
			my $contigName = shift @temp_arr;
		    my $contigSeq = join("", @rows);
		    $contigSeq =~ s/\s//g; #remove white space
		    $seq_hash_ref->{$contigName} = $contigSeq;
			$seq_desc_hash_ref->{$contigName} = join(" ", @temp_arr);
#			print " name = \\$contigName\\, seq  = \\$contigSeq\\\n\n";
		}
    }

    #reset the read seperator
    $/ = $oldseperator;
}

