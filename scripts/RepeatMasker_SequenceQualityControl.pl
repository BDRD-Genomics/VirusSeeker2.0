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

my $usage = '
This script will check each .masked file in the given directory.
Some sequences have only/lots of Ns because masked by RepeatMasker.
1) Sequences that do not have greater than 50 nt of consecutive
sequence without N will be put into file .cdhit_out.masked.badSeq 
2) Sequences with >= 40% of total length of being masked will be put 
into file .cdhit_out.masked.RepeatLowComplexSeq
3) output original sequnce into .cdhit_out.masked.goodSeq file.

perl script <sample dir>
<sample dir> = full path of the folder holding files for this sequence run
               without last "/"

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;
my $percent_masked_cutoff = 0.4;

####################################################################################################
# get name of the sample and path to the data
if ($dir =~/(.+)\/$/) {
	$dir = $1;
}
my $sample_name = (split(/\//,$dir))[-1];
print $sample_name,"\n";

my $original_fa_file = $dir."/".$sample_name.".allSeqs.mmseqs_rep_seq.fasta";
my %original_fa_seq = ();
&read_FASTA_data($original_fa_file, \%original_fa_seq);

my $total_seq = 0;
my $good_seq = 0;
my $bad_seq = 0;
my $RepeatLowComplex_seq = 0;
my $OutFile1 = $dir."/".$sample_name.".goodSeq.fa";
my $OutFile2 = $dir."/".$sample_name.".badSeq.fa";
my $OutFile3 = $dir."/".$sample_name.".RepeatLowComplexSeq.fa";

open (OUT1, ">$OutFile1") or die "can not open $OutFile1\n";
open (OUT2, ">$OutFile2") or die "can not open $OutFile2\n";
open (OUT3, ">$OutFile3") or die "can not open $OutFile3\n";

my $full_path = $dir."/".$sample_name."_RepeatMasker"; # RepeatMasker directory
opendir(SubDH, $full_path) or die "can not open dir $full_path!\n";
foreach my $file (readdir SubDH) {

	if ($file =~ /\.masked$/) { # masked sequence
		my $maskedFile = $full_path."/".$file;
		my %seq = ();
		&read_FASTA_data($maskedFile, \%seq);

		# check for contiguous bases >= 50 bp (non-Ns) 
		foreach my $read_id (keys %seq) {
			$total_seq++;
			my $seq_temp = $seq{$read_id};
			my $goodQuality=$seq_temp=~/[ACTG]{50,}/; 
			if($goodQuality) {
				my $length_masked = ($seq_temp =~ tr/N/N/);
				my $length_total = length $seq_temp;
                my $percent_masked = $length_masked/$length_total;
#				print  ">$read_id\n";
#				print  $seq{$read_id}, "\n";
#	            print "total length $length_total, total number of Ns $length_masked, percentage $percent_masked\n";
				if ($percent_masked >= $percent_masked_cutoff) { # output low complex sequence
					print OUT3 ">$read_id\n";
					print OUT3 $seq{$read_id}, "\n";
					$RepeatLowComplex_seq++;
				}
				else { # output good seq
					print OUT1 ">$read_id\n";
					print OUT1 $seq{$read_id}, "\n";
#					print OUT1 $original_fa_seq{$read_id}, "\n";

					$good_seq++;
				}
			}
			else { # output bad sequence
				print OUT2 ">$read_id\n";
				print  OUT2 "$seq{$read_id}\n";
				$bad_seq++;
			}
		}
	}
}

print OUT2 "total seq = $total_seq\n";
print OUT2 "good seq = $good_seq\n";
print OUT2 "bad seq = $bad_seq\n";
print OUT2 "Repeat and Low complexicity seq = $RepeatLowComplex_seq\n";


print OUT3 "total seq = $total_seq\n";
print OUT3 "Repeat and Low complexicity seq = $RepeatLowComplex_seq\n";

close(OUT1);
close(OUT2);
close(OUT3);

exit;

############################################################################
sub read_FASTA_data () {
    my ($fastaFile, $hash_ref) = @_;

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
		    my $seqName = shift @rows;
			my @temp = split (/\s/, $seqName);
			$seqName = shift @temp;
		    my $Seq = join("", @rows);
		    $Seq =~ s/\s//g; #remove white space
		    $hash_ref->{$seqName} = $Seq;
#			print "name = $seqName\n";
#			print "seq = \\$Seq\\\n";
		}
    }

	close FastaFile;
    #reset the read seperator
    $/ = $oldseperator;
}
