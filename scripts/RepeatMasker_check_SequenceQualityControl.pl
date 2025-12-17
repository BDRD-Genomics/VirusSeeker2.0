###############################################################################
# Author: BDRD <usn.detrick.nmrc.mbx.genomics-reach-back@health.mil>
#
# License:
# VirusSeeker 2.0 is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3 of the License or any later version. For 
# details, please refer to https://www.gnu.org/licenses/
###############################################################################


#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample dir>
<sample dir> = full path of the folder holding files for a sample

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

if ($dir =~/(.+)\/$/) {
	$dir = $1;
}
my $sample_name = (split(/\//,$dir))[-1];
#print $sample_name,"\n";

my $finished = &check_QC_read_number($dir);
#print $finished; 
exit ($finished);

##########################################################################
sub check_QC_read_number {
	my ( $dir ) = @_;
	my $tot_cdhit_seq=0;
	my $tot_seq = 0;

	my $in_file = $dir."/".$sample_name.".allSeqs.mmseqs_rep_seq.fasta";
	print "in file is $in_file\n";
	$tot_cdhit_seq = `grep \">\" $in_file |wc -l`;

	my $badSeq_file = $dir."/".$sample_name.".badSeq.fa";
	print "QC output file is $badSeq_file\n";
	my $temp = `grep \"total seq\" $badSeq_file`;
	if ($temp =~ /total seq = (\d+)/) {
		$tot_seq = $1;
	}

	print "total unique seq in input file: $tot_cdhit_seq\n";
	print "total unique sequence in QC output file: $tot_seq\n\n"; 

	if($tot_cdhit_seq==$tot_seq) { return 0; } 
	else { 
		return 10; 
	}
}
