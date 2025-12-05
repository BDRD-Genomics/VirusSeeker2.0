###################################################################
# Author: Kyle Long <kyle.a.long8.ctr@mail.mil>
# Additional Contact: Regina Cer <regina.z.cer.civ@mail.mil>
#
# License:
# VirusSeeker 2.0 is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3 of the License or any later version. For 
# details, please refer to https://www.gnu.org/licenses/
###############################################################################


#!/usr/bin/perl -w
use strict;
my $usage="
perl $0 <sample dir>
<sample dir> = full path of the folder holding files for this sample
               without last \"/\"
";

die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

if ($dir =~/(.+)\/$/) {
        $dir = $1;
}
my $sample_name = (split(/\//,$dir))[-1];
#print $sample_name,"\n";

my $finished = &check_split_output($dir);
#print $finished; 
exit ($finished);

##############################################################
sub check_split_output {
	my ( $dir ) = @_;
	my $tot_in_seq=0;
	my $tot_out_seq = 0;

	my $in_file = $dir."/".$sample_name.".RefGenomeFiltered.goodSeq_and_top3.fastq";
	my $total_lines = `wc -l < $in_file`;
	$tot_in_seq = $total_lines / 4;

	my $Assembly_dir_SE = $dir."/".$sample_name."_Assembly_SE";
	opendir(DH, $Assembly_dir_SE) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fastq$/) { 
			my $in = $Assembly_dir_SE."/".$file;
			my $total_lines = `wc -l < $in`;
			my $temp = $total_lines / 4;
			$tot_out_seq += $temp;                    
		}
	}
	close DH; 

	print "total input sequence: $tot_in_seq\n";
	print "total output sequence : $tot_out_seq\n"; 

	if($tot_in_seq==$tot_out_seq) { return 0; }  # use 0 to be consistent with Linux convention
	else { # remove content of assembly directory, to rerun split 
		opendir(DH, $Assembly_dir_SE) or return 10;
		foreach my $file (readdir DH) {	
			my $faFile = $Assembly_dir_SE."/".$file;
#            print "$faFile\n"; 
			unlink $faFile; 
		}
	  	close DH;  
	  	return 10; 
	}
}
