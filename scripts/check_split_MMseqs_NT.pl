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

my $finished = &check_split_MegaBLAST($dir);
#print $finished; 
exit ($finished);

##########################################################################
sub check_split_MegaBLAST {
	my ( $dir ) = @_;
	my $tot_in_seq=0;
	my $tot_out_seq = 0;

	my $in_file = $dir."/".$sample_name.".VIRUSDB_HIT.fa";
	print "in file is $in_file\n";
	$tot_in_seq = `grep \">\" $in_file |wc -l`;

	my $MegaBLAST_dir = $dir."/".$sample_name."_MMseqs_NT";
	opendir(DH, $MegaBLAST_dir) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fasta$/ ) {
			my $faFile = $MegaBLAST_dir."/".$file;
			my $temp = `grep \">\" $faFile |wc -l`;
			$tot_out_seq += $temp;
		}
	}
	close DH;

	print "total seq in input file: $tot_in_seq\n";
	print "total sequence in splited output dir: $tot_out_seq\n\n"; 

	if($tot_in_seq==$tot_out_seq) { return 0; } # use 0 to be consistant with Linux convention. 
	else { 
		return 10; 
	}
}
