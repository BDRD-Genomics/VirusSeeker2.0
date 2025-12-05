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


#!/usr/bin/perl
use strict;

my $usage = '
This script accepts one contig fasta file and one singleton fasta file 
and put all the sequences into one file. The combined data set naming 
convention is: sequences from contig file named contig1, contig2 ...
This ensures each sequence in the output file has a unique name.
Sequences from singleton file already has unique names. 

perl script <input contig fasta file><input singleton fasta file>
<input contig fasta file> = sequence file with full path 

output to standard output.

';
die $usage unless scalar @ARGV == 1;
my ($F_contig) = @ARGV;

# For the contig sequence file
my $i = 1;
#keep old read seperator and set new read seperator to ">"
my $oldseperator = $/;
$/ = ">";
open (IN, $F_contig) or die "Can't Open FASTA file: $F_contig!\n";
while (my $line = <IN>){
	chomp $line;

	# Discard blank lines
	if ($line =~ /^\s*$/) {
		next;
	}	
	# discard the first line which only has ">", keep the rest
	elsif ($line ne ">") {
		print ">FirstAssembly_contig", $i, " ", $line;
		$i++ 
	}
}

#reset the read seperator
$/ = $oldseperator;
close IN; 

exit;

