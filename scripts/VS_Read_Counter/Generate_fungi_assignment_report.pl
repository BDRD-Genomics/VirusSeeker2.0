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

perl script <sample dir> 
<sample dir> = full path to the directory holding files for the given 
               library without the last "/"
               e.g. ~/tools/Illumina_VirusDiscoveryPipeline/data/MiSeq_run_1/I10_12310_Project404_CSF_Glaser_Encephalitis_V11T00919_TruSeq-AD007_CAGATC

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

my @temp = split("\/", $dir);
my $sample_name = pop @temp;
# print "lib is $sample_name\n";

####################################
# read in original sequences
my $fasta_file = $dir."/".$sample_name.".allSeqs.mmseqs_all_seqs.fasta"; #

my %seq = (); # read_ID => sequence
my %seq_desc = (); # read_ID => description
&read_FASTA_data($fasta_file, \%seq, \%seq_desc);

####################################
# output files 
my $out3 = $dir."/".$sample_name.".FungiReadsAssignmentReport";
open (OUT3, ">$out3") or die "can not open file $out3!\n";
my $out4 = $dir."/".$sample_name.".FungiReads_all.fa";
open (OUT4, ">$out4") or die "can not open file $out4!\n";

##############################################################################
# variable to keep all the blast output files
my @blast_files_MegaBlast_NT = (); # all MegaBlast output files
my @blast_files_blastn = (); # all blastn.out files
my @blast_files_blastx = (); # all blastx.out files

# all the viral read ID
my %viral_reads_MegaBlast_NT = ();
my %viral_reads_blastn = ();
my %viral_reads_blastx = ();
my %viral_reads_blastx_VIRUSDB = ();

# viral_lineage => number of reads assigned to this lineage in the sample
my %viral_lineage_to_num_reads = ();
my %viral_reads_blast_info = ();    # readID => information about this read
my %viral_reads_in_lineage_MegaBlast_NT = ();    # lineage => [read ID]
my %viral_reads_in_lineage_blastn = ();    # lineage => [read ID]
my %viral_reads_in_lineage_blastx = ();    # lineage => [read ID]

# all the fungi viral read ID
my %fungi_reads_MegaBlast_NT = ();
my %fungi_reads_blastn = ();
my %fungi_reads_blastx = ();

# fungi reads viral lineage => number of reads assigned to this lineage in the sample
my %fungi_reads_lineage_to_num_reads = ();
my %fungi_reads_blast_info = ();    # readID => information about this read
my %fungi_reads_in_lineage_MegaBlast_NT = ();    # lineage => [read ID]
my %fungi_reads_in_lineage_blastn = ();    # lineage => [read ID]
my %fungi_reads_in_lineage_blastx = ();    # lineage => [read ID]

my @unassigned_reads = ();
my @phage_reads = ();
##################################
# generate unassigned reads file
my $BLASTX_NR_dir = $dir."/".$sample_name."_BLASTX_NR";



# category => num of sequence assigned to this category by blastn
my %Assignment_blastn = ( 
	"Fungi" => 0,
	"Fungi" => 0,
	"Homo" => 0,
	"Mus" => 0,
	"Phage" => 0,
	"Viruses" => 0,
	"other" => 0,
	"unassigned" => 0,
	"Fungi" => 0,
);

# category => num of sequence assigned to this category by MegaBLAST 
my %Assignment_MegaBlast_NT = ();
foreach my $key (keys %Assignment_blastn) {
	$Assignment_MegaBlast_NT{$key} = 0;
}
# category => num of sequence assigned to this category by blastx 
my %Assignment_blastx = ();
foreach my $key (keys %Assignment_blastn) {
	$Assignment_blastx{$key} = 0;
}

# category => num of sequence assigned to this category by blastx VIRUSDB 
my %Assignment_blastx_VIRUSDB = ();
foreach my $key (keys %Assignment_blastn) {
	$Assignment_blastx_VIRUSDB{$key} = 0;
}


####################################################################################
# to obtain viral read information
# MegaBLAST NT step
my $MegaBLAST_dir = $dir."/".$sample_name."_MegaBLAST_NT";
# enter directory where blast results resides
opendir (DIR, $MegaBLAST_dir) or die "can not open dir $MegaBLAST_dir!\n";
foreach my $blast_file (readdir DIR) {
	if ($blast_file =~ /megablast\.parsed$/) {
		# print "blastn parsed file $blast_file\n";
		my $parsed = $MegaBLAST_dir."/".$blast_file;
		my $blast_out = $blast_file;
		$blast_out =~ s/\.megablast\.parsed/\.megablast\.out/;
		$blast_out = $MegaBLAST_dir."/".$blast_out;
		push @blast_files_MegaBlast_NT, $blast_out;

		&collect_information($parsed, \%Assignment_MegaBlast_NT, \%viral_reads_MegaBlast_NT, \%viral_reads_in_lineage_MegaBlast_NT, \%viral_lineage_to_num_reads, \%fungi_reads_MegaBlast_NT, \%fungi_reads_in_lineage_MegaBlast_NT, \%fungi_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%fungi_reads_blast_info, \@phage_reads );
	}
}
closedir DIR;


# BLASTN step
my $BLASTN_dir = $dir."/".$sample_name."_BLASTN_NT";
# enter directory where blastn results resides
opendir (BNDIR, $BLASTN_dir) or die "can not open dir $BLASTN_dir!\n";
foreach my $blast_file (readdir BNDIR) {
	if ($blast_file =~ /blastn\.parsed$/) {
		# print "blastn parsed file $blast_file\n";
		my $blast_out = $blast_file;
		$blast_out =~ s/\.blastn\.parsed/\.blastn\.out/;
		$blast_out = $BLASTN_dir."/".$blast_out;
		push @blast_files_blastn, $blast_out;
		my $parsed = $BLASTN_dir."/".$blast_file;
		&collect_information($parsed, \%Assignment_blastn, \%viral_reads_blastn, \%viral_reads_in_lineage_blastn, \%viral_lineage_to_num_reads, \%fungi_reads_blastn, \%fungi_reads_in_lineage_blastn, \%fungi_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%fungi_reads_blast_info );
	}
}
closedir BNDIR;

# BLASTX NR step
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
		&collect_information($parsed, \%Assignment_blastx, \%viral_reads_blastx, \%viral_reads_in_lineage_blastx, \%viral_lineage_to_num_reads, \%fungi_reads_blastx, \%fungi_reads_in_lineage_blastx, \%fungi_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%fungi_reads_blast_info );
	}
}
closedir BXDIR;


my $c1 = "############################################################\n\n";


# generate fungi reads report
print OUT3 $c1;
foreach my $lineage (sort {$fungi_reads_lineage_to_num_reads{$a} <=> $fungi_reads_lineage_to_num_reads{$b}} keys %fungi_reads_lineage_to_num_reads) {
	print OUT3 $lineage, "\ttotal number of reads: ", $fungi_reads_lineage_to_num_reads{$lineage}, "\n\n";
	print OUT3 "QueryName\tQuerylength\t         HitName       \tHitLen\t                             HitDesc                       \tAlnLen\t%ID\tHitStart\tHitEnd\te\n";

	if (defined $fungi_reads_in_lineage_MegaBlast_NT{$lineage}) {
		if (scalar @{$fungi_reads_in_lineage_MegaBlast_NT{$lineage}}) {
			print OUT3 "reads from MegaBlast NT:\n";
			foreach my $read (sort {$a cmp $b} @{$fungi_reads_in_lineage_MegaBlast_NT{$lineage}}) {
				print OUT3 $fungi_reads_blast_info{$read}, "\n";
			}
		}
	}


	if (defined $fungi_reads_in_lineage_blastn{$lineage}) {
		if (scalar @{$fungi_reads_in_lineage_blastn{$lineage}}) {
			print OUT3 "reads from blastn:\n";
			foreach my $read (sort {$a cmp $b} @{$fungi_reads_in_lineage_blastn{$lineage}}) {
				print OUT3 $fungi_reads_blast_info{$read}, "\n";
			}
		}
		print OUT3 "\n";
	}

	if (defined $fungi_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$fungi_reads_in_lineage_blastx{$lineage}}) {
			print OUT3 "reads from blastx:\n";
			foreach my $read (sort {$a cmp $b} @{$fungi_reads_in_lineage_blastx{$lineage}}) {
				print OUT3 $fungi_reads_blast_info{$read}, "\n";
			}
		}
		print OUT3 "\n";
	}


	##########################################################################
	# print out sequences
	if (defined $fungi_reads_in_lineage_MegaBlast_NT{$lineage}) {
		if (scalar @{$fungi_reads_in_lineage_MegaBlast_NT{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$fungi_reads_in_lineage_MegaBlast_NT{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

	if (defined $fungi_reads_in_lineage_blastn{$lineage}) {
		if (scalar @{$fungi_reads_in_lineage_blastn{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$fungi_reads_in_lineage_blastn{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

	if (defined $fungi_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$fungi_reads_in_lineage_blastx{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$fungi_reads_in_lineage_blastx{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

	print OUT3 $c1;
}


###############################################################################
# get all the fungi reads and put the sequence into output file:
foreach my $lineage (keys %fungi_reads_lineage_to_num_reads) {
	foreach my $read (@{$fungi_reads_in_lineage_MegaBlast_NT{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}

	foreach my $read (@{$fungi_reads_in_lineage_blastn{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}
	foreach my $read (@{$fungi_reads_in_lineage_blastx{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}
#	foreach my $read (@{$fungi_reads_in_lineage_blastx_VIRUSDB{$lineage}}) {
#		print OUT4 ">$read $seq_desc{$read}\n";
#		print OUT4 $seq{$read}, "\n";
#	}
}

print OUT1 "# Finished Assignment Report\n";

exit;

#####################################################################################
# collecte information from given BLAST parsed output file
sub collect_information {
	################################## modified 04/2014
	my ($blast_parsed_file, $category_hash_ref, $viral_reads_hash_ref, $viral_reads_in_lineage_hash_ref, $viral_lineage_to_num_reads_hash_ref,  $fungi_reads_hash_ref, $fungi_reads_in_lineage_hash_ref, $fungi_reads_lineage_to_num_reads_hash_ref, $unassigned_reads_arr_ref, $viral_reads_blast_info_ref, $fungi_reads_blast_info_ref ) = @_;

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
			case "Fungi" { $category_hash_ref->{"Fungi"}++	}
			case "Fungi" { $category_hash_ref->{"Fungi"}++ }
			case "Homo" { $category_hash_ref->{"Homo"}++ }
			case "Mus" { $category_hash_ref->{"Mus"}++ }
			case "Phage" {$category_hash_ref->{"Phage"}++ }
			case "Viruses" { $category_hash_ref->{"Viruses"}++ }
			case "other" {$category_hash_ref->{"other"}++ }
			case "unassigned" {$category_hash_ref->{"unassigned"}++}
			case "Fungi" {$category_hash_ref->{"Fungi"}++ } 
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
		}elsif ($category eq "Fungi"){
			$fungi_reads_hash_ref->{$read_ID} = 1;
			$fungi_reads_blast_info_ref->{$read_ID} = $desc;
			if (!(defined $fungi_reads_in_lineage_hash_ref->{$lineage})) {
				$fungi_reads_in_lineage_hash_ref->{$lineage} = [$read_ID];
			}
			else {
				push @{$fungi_reads_in_lineage_hash_ref->{$lineage}}, $read_ID;
			}
		
			if (defined $fungi_reads_lineage_to_num_reads_hash_ref->{$lineage}) {
				$fungi_reads_lineage_to_num_reads_hash_ref->{$lineage}++;
			}
			else {
				$fungi_reads_lineage_to_num_reads_hash_ref->{$lineage} = 1;
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

    # check
#	 foreach my $key (keys %fastaSeq){
#	      print "Here is the key for fasta seq: $key \t $fastaSeq{$key}\n";
#	 }

    #reset the read seperator
    $/ = $oldseperator;
}

