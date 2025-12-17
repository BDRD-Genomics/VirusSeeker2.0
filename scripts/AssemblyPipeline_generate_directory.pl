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

my $usage = "
This script will create the appropriate directory for the samples in a given directory.

perl $0  <input run folder><output assembly directory> 
<run folder> = full path of the folder holding files to be analysized
               without the last \"/\"
<output assembly directory> = path to the assembly directory

";

die $usage unless scalar @ARGV == 2;
my ( $in_dir , $out_dir ) = @ARGV;

opendir(DH, $in_dir) or die "Can not open dir $in_dir!\n";
foreach my $name (readdir DH) {
	if (!($name eq ".") && !($name eq "..")) {
		my $full_path = $in_dir."/".$name;
		if (-d $full_path) { # is a directory
			my $com = "cp --parents $full_path/*fastq.gz $out_dir/ \n";
			print $com, "\n";
			system( $com );

			my $com = "cp --parents $full_path/*RemoveAdapter.fastq $out_dir/ \n";
			print $com, "\n";
			system( $com );

			my $com = "cp --parents $full_path/*RemoveAdapter.report $out_dir/ \n";
			print $com, "\n";
			system( $com );

			my $com = "cp --parents $full_path/*_adapter.txt $out_dir/ \n";
			print $com, "\n";
			system( $com );


			my $com = "cp --parents $full_path/status_log/j1_RemoveAdapter* $out_dir/ \n";
			print $com, "\n";
			system( $com );
		}
	}
}

exit;

