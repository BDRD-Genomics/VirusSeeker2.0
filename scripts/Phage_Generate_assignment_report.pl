###############################################################################
# Author: Kyle Long <kyle.a.long8.ctr@health.mil>
# Additional Contact: Regina Cer <regina.z.cer.civ@health.mil>
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
generate a report for all the phage sequences. 

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
my $fasta_file = $dir."/".$sample_name.".allSeqs.mmseqs_rep_seq.fasta";
my %seq = (); # read_ID => sequence
my %seq_desc = (); # read_ID => description
&read_FASTA_data($fasta_file, \%seq, \%seq_desc);

####################################
# output files 
my $out1 = $dir."/".$sample_name.".PhageAssignmentReport";
open (OUT1, ">$out1") or die "can not open file $out1!\n";

my $out2 = $dir."/".$sample_name.".phiX174.fasta";
open (OUT2, ">$out2") or die "can not open file $out2!\n";

my $out3 = $dir."/".$sample_name.".PhageReads_All.fa";
open (OUT3, ">$out3") or die "can not open file $out3!\n";


##############################################################################
# variable to keep all the blast output files
my @blast_files_MMseqs_VIRUSDB = (); # all MMseqs_VIRUSDB output files
my @blast_files_BLASTx_VIRUSDB = (); # all BLASTx_VIRUSDB output files
my @blast_files_MMseqs_NT = (); # all MMseqs output files
#my @blast_files_blastn = (); # all blastn.out files
my @blast_files_blastx = (); # all blastx.out files

# all the phage read ID: significant hits at BLASTn VIRUSDB, BLASTx VIRUSDB,
# MMseeqs_NT, and blastx NR 
my %phage_reads_MMseqs_VIRUSDB = ();
my %phage_reads_BLASTx_VIRUSDB = ();
my %phage_reads_MMseqs_NT = ();
#my %phage_reads_blastn = ();
my %phage_reads_blastx = ();

# phage_lineage => number of reads assigned to this lineage in the sample
my %phage_lineage_to_num_reads = ();
my %phage_reads_blast_info =();    # readID => information about this read
my %phage_reads_in_lineage_MMseqs_VIRUSDB = ();    # lineage => [read ID]
my %phage_reads_in_lineage_BLASTx_VIRUSDB = ();    # lineage => [read ID]
my %phage_reads_in_lineage_MMseqs_NT = ();    # lineage => [read ID]
#my %phage_reads_in_lineage_blastn = ();    # lineage => [read ID]
my %phage_reads_in_lineage_blastx = ();    # lineage => [read ID]

####################################################################################
# to obtain phage read information
# BLASTn VIRUSDB step
my $MMseqs_VIRUSDB_dir = $dir."/".$sample_name."_MMseqs_VIRUSDB";
# enter directory where blastn VIRUSDB results resides
opendir (BNDIR, $MMseqs_VIRUSDB_dir) or die "can not open dir $MMseqs_VIRUSDB_dir!\n";
foreach my $blast_file (readdir BNDIR) {
	if ($blast_file =~ /mmseqsv\.parsed$/) {
		#print "blastn VIRUSDB parsed file $blast_file\n";
		my $blast_out = $blast_file;
		$blast_out =~ s/\.mmseqsv\.parsed/\.mmseqsv\.out/;
		$blast_out = $MMseqs_VIRUSDB_dir."/".$blast_out;
		push @blast_files_MMseqs_VIRUSDB, $blast_out;
		my $parsed = $MMseqs_VIRUSDB_dir."/".$blast_file;
		&collect_phage_information($parsed, \%phage_reads_MMseqs_VIRUSDB, \%phage_reads_in_lineage_MMseqs_NT, \%phage_lineage_to_num_reads, \%phage_reads_blast_info );
	}
}
closedir BNDIR;

# BLASTX VIRUSDB step
my $BLASTX_VIRUSDB_dir = $dir."/".$sample_name."_BLASTX_VIRUSDB";
# enter directory where blastx VIRUSDB results resides
opendir (BXDIR, $BLASTX_VIRUSDB_dir) or die "can not open dir $BLASTX_VIRUSDB_dir!\n";
foreach my $blast_file (readdir BXDIR) {
	if ($blast_file =~ /BLASTX_VIRUSDB\.parsed$/) {
		# print "blast parsed file $blast_file\n";
		my $blast_out = $blast_file;
		$blast_out =~ s/BLASTX_VIRUSDB\.parsed/BLASTX_VIRUSDB\.out/;
		$blast_out = $BLASTX_VIRUSDB_dir."/".$blast_out;
		push @blast_files_blastx, $blast_out;
		my $parsed = $BLASTX_VIRUSDB_dir."/".$blast_file;
		&collect_phage_information($parsed, \%phage_reads_BLASTx_VIRUSDB, \%phage_reads_in_lineage_blastx, \%phage_lineage_to_num_reads, \%phage_reads_blast_info);
	}
}
closedir BXDIR;


# MMseqs NT step
my $MMseqs_dir = $dir."/".$sample_name."_MMseqs_NT";
# enter directory where blast results resides
if ( -e $MMseqs_dir) { # some times this does not exist
	opendir (DIR, $MMseqs_dir) or die "can not open dir $MMseqs_dir!\n";
	foreach my $blast_file (readdir DIR) {
		if ($blast_file =~ /mmseqs\.parsed$/) {
			# print "blastn parsed file $blast_file\n";
			my $parsed = $MMseqs_dir."/".$blast_file;
			my $blast_out = $blast_file;
			$blast_out =~ s/\.mmseqs\.parsed/\.mmseqs\.out/;
			$blast_out = $MMseqs_dir."/".$blast_out;
			push @blast_files_MMseqs_NT, $blast_out;

			&collect_phage_information($parsed, \%phage_reads_MMseqs_NT, \%phage_reads_in_lineage_MMseqs_NT, \%phage_lineage_to_num_reads, \%phage_reads_blast_info );
		}
	}
}
closedir DIR;

# # BLASTN step
# my $BLASTN_dir = $dir."/".$sample_name."_BLASTN_NT";
# # enter directory where blastn results resides
# if ( -e $BLASTN_dir) {
# 	opendir (BNDIR, $BLASTN_dir) or die "can not open dir $BLASTN_dir!\n";
# 	foreach my $blast_file (readdir BNDIR) {
# 		if ($blast_file =~ /blastn\.parsed$/) {
# 			# print "blastn parsed file $blast_file\n";
# 			my $blast_out = $blast_file;
# 			$blast_out =~ s/\.blastn\.parsed/\.blastn\.out/;
# 			$blast_out = $BLASTN_dir."/".$blast_out;
# 			push @blast_files_blastn, $blast_out;
# 			my $parsed = $BLASTN_dir."/".$blast_file;
# 			&collect_phage_information($parsed, \%phage_reads_blastn, \%phage_reads_in_lineage_blastn, \%phage_lineage_to_num_reads, \%phage_reads_blast_info);
# 		}
# 	}
# }
# closedir BNDIR;

# BLASTX NR step
my $BLASTX_dir = $dir."/".$sample_name."_BLASTX_NR";
# enter directory where blastx results resides
if ( -e $BLASTX_dir) {
	opendir (BXDIR, $BLASTX_dir) or die "can not open dir $BLASTX_dir!\n";
	foreach my $blast_file (readdir BXDIR) {
		if ($blast_file =~ /blastx\.parsed$/) {
			# print "blast parsed file $blast_file\n";
			my $blast_out = $blast_file;
			$blast_out =~ s/\.blastx\.parsed/\.blastx\.out/;
			$blast_out = $BLASTX_dir."/".$blast_out;
			push @blast_files_blastx, $blast_out;
			my $parsed = $BLASTX_dir."/".$blast_file;
			&collect_phage_information($parsed, \%phage_reads_blastx, \%phage_reads_in_lineage_blastx, \%phage_lineage_to_num_reads, \%phage_reads_blast_info);
		}
	}
}
closedir BXDIR;

#####################################################################################
# print out report for this library
print OUT1 $dir, "\n";

##############################################################################
# phage lineage, number of phage reads in this lineage, BLAST alginment statistics.

my $c1 = "############################################################\n\n";
print OUT1 $c1;

foreach my $lineage (sort {$phage_lineage_to_num_reads{$a} <=> $phage_lineage_to_num_reads{$b}} keys %phage_lineage_to_num_reads) {
	# print out phiX174 sequences
	if ($lineage =~ /Enterobacteria phage phiX174 sensu lato/) {
		# print out phiX174 sequences
		if (defined $phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}}) {
					print OUT2 ">$read $seq_desc{$read}\n";
					print OUT2 $seq{$read}, "\n";
				}
			}
		}

		if (defined $phage_reads_in_lineage_BLASTx_VIRUSDB{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_blastx{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastx{$lineage}}) {
					print OUT2 ">$read $seq_desc{$read}\n";
					print OUT2 $seq{$read}, "\n";
				}
			}
		}

		if (defined $phage_reads_in_lineage_MMseqs_NT{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_MMseqs_NT{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_MMseqs_NT{$lineage}}) {
					print OUT2 ">$read $seq_desc{$read}\n";
					print OUT2 $seq{$read}, "\n";
				}
			}
		}

		# if (defined $phage_reads_in_lineage_blastn{$lineage}) {
		# 	if (scalar @{$phage_reads_in_lineage_blastn{$lineage}}) {
		# 		foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastn{$lineage}}) {
		# 			print OUT2 ">$read $seq_desc{$read}\n";
		# 			print OUT2 $seq{$read}, "\n";
		# 		}
		# 	}
		# }

		if (defined $phage_reads_in_lineage_blastx{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_blastx{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastx{$lineage}}) {
					print OUT2 ">$read $seq_desc{$read}\n";
					print OUT2 $seq{$read}, "\n";
				}
			}
		}
	}
	else { # for other lineages
		print OUT1 $lineage, "\ttotal number of reads: ", $phage_lineage_to_num_reads{$lineage}, "\n\n";
		print OUT1 "QueryName\tQuerylength\tHitName\tHitLen\tHitDesc\tAlnLen\t%ID\tHitStart\tHitEnd\te\n";
		if (defined $phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}}) {
				print OUT1 "reads from MMseqs VIRUSDB:\n";
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}}) {
					print OUT1 $phage_reads_blast_info{$read}, "\n";
				}
			}
			print OUT1 "\n";
		}

		if (defined $phage_reads_in_lineage_BLASTx_VIRUSDB{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_BLASTx_VIRUSDB{$lineage}}) {
				print OUT1 "reads from BLASTx VIRUSDB:\n";
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_BLASTx_VIRUSDB{$lineage}}) {
					print OUT1 $phage_reads_blast_info{$read}, "\n";
				}
			}
			print OUT1 "\n";
		}

		if (defined $phage_reads_in_lineage_MMseqs_NT{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_MMseqs_NT{$lineage}}) {
				print OUT1 "reads from MMSeqs NT:\n";
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_MMseqs_NT{$lineage}}) {
					print OUT1 $phage_reads_blast_info{$read}, "\n";
				}
			}
			print OUT1 "\n";
		}

		# if (defined $phage_reads_in_lineage_blastn{$lineage}) {
		# 	if (scalar @{$phage_reads_in_lineage_blastn{$lineage}}) {
		# 		print OUT1 "reads from BLASTn NT:\n";
		# 		foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastn{$lineage}}) {
		# 			print OUT1 $phage_reads_blast_info{$read}, "\n";
		# 		}
		# 	}
		# 	print OUT1 "\n";
		# }

		if (defined $phage_reads_in_lineage_blastx{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_blastx{$lineage}}) {
				print OUT1 "reads from BLASTx NR:\n";
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastx{$lineage}}) {
					print OUT1 $phage_reads_blast_info{$read}, "\n";
				}
			}
			print OUT1 "\n";
		}

		# print out phage sequences in each taxonomy lineage
		if (defined $phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_MMseqs_VIRUSDB{$lineage}}) {
					print OUT3 ">$read $seq_desc{$read}\n";
					print OUT3 $seq{$read}, "\n";

				}
			}
		}

		if (defined $phage_reads_in_lineage_BLASTx_VIRUSDB{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_blastx{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastx{$lineage}}) {
					print OUT3 ">$read $seq_desc{$read}\n";
					print OUT3 $seq{$read}, "\n";
				}
			}
		}

		if (defined $phage_reads_in_lineage_MMseqs_NT{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_MMseqs_NT{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_MMseqs_NT{$lineage}}) {
					print OUT3 ">$read $seq_desc{$read}\n";
					print OUT3 $seq{$read}, "\n";
				}
			}
		}

		# if (defined $phage_reads_in_lineage_blastn{$lineage}) {
		# 	if (scalar @{$phage_reads_in_lineage_blastn{$lineage}}) {
		# 		foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastn{$lineage}}) {
		# 			print OUT3 ">$read $seq_desc{$read}\n";
		# 			print OUT3 $seq{$read}, "\n";
		# 		}
		# 	}
		# }

		if (defined $phage_reads_in_lineage_blastx{$lineage}) {
			if (scalar @{$phage_reads_in_lineage_blastx{$lineage}}) {
				foreach my $read (sort {$a cmp $b} @{$phage_reads_in_lineage_blastx{$lineage}}) {
					print OUT3 ">$read $seq_desc{$read}\n";
					print OUT3 $seq{$read}, "\n";
				}
			}
		}

		print OUT1 $c1;
	}
}


print OUT1 "# Finished Assignment Report\n";

close OUT1;
close OUT2;
close OUT3;
exit;

#####################################################################################
# collecte information from given BLAST parsed output file
sub collect_phage_information {
	my ($blast_parsed_file, $phage_reads_hash_ref, $phage_reads_in_lineage_hash_ref, $phage_lineage_to_num_reads_hash_ref, $phage_reads_blast_info_ref) = @_;

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

#		print "readID = \\$read_ID\\, length = \\$length\\, category = \\$category\\, lineage = \\$lineage\\, desc = \\$desc\\\n\n";
		$desc = $read_ID."\t".$length."\t".$desc; 
		if ($category eq "Phage") {
			$phage_reads_hash_ref->{$read_ID} = 1;
			$phage_reads_blast_info_ref->{$read_ID} = $desc;
			if (!(defined $phage_reads_in_lineage_hash_ref->{$lineage})) {
				$phage_reads_in_lineage_hash_ref->{$lineage} = [$read_ID];
			}
			else {
				push @{$phage_reads_in_lineage_hash_ref->{$lineage}}, $read_ID;
			}
		
			if (defined $phage_lineage_to_num_reads_hash_ref->{$lineage}) {
				$phage_lineage_to_num_reads_hash_ref->{$lineage}++;
			}
			else {
				$phage_lineage_to_num_reads_hash_ref->{$lineage} = 1;
			}
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

