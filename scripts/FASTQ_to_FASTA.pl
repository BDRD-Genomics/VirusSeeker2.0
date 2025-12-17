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


# in this version the read name ID is trimmed for "/2" for the purpose 
# to fit the pipeline because I used bamToFastq to generate the 
# fastq from bam mapping file and the read ID were postfixed with "/2" 


#!/usr/bin/perl -w

use strict;

my $usage = '
Given a fastq file, this script will conver it to a fasta format file.

perl script <fastq file>
<fastq file> = full path to the fastq file

output to standard output.

';
die $usage unless scalar @ARGV == 1;
my ($inFile) = @ARGV;

open (IN, "<$inFile") or die "Can't Open file: $inFile";

while (my $line = <IN>){
	my $info = $line; # first line
	chomp $line;
	my $name = $line;
	$name =~ s/@//;
#	chop $name; #commented out by CWP
#	chop $name; #commented out by CWP
	print ">$name\n";
	$line = <IN>; # second line: sequence
	print "$line";
	$line = <IN>; # third line: + 
	$line = <IN>; # fourth line: quality line
}

close IN;

exit;

