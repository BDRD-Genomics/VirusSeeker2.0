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



#!/usr/bin/perl -w

use strict;

my $usage = '
Given the fastq-join output un1, un2 and join fastq files, the script will
rename the id in the sequencing data. Joined reads will have a prefix of
"join", un1 reads (unjoined SE1 reads) will have a prefix of "un1" and un2 
reads will have a prefix of "un2".

perl script < full path and the fastq-join output file prefix for un1, un2 
and join fastq files>

';

die $usage unless scalar @ARGV == 1;
my ($inName) = @ARGV;

#my $file_join=$inName.".join.fq";
#my $file_un1=$inName.".un1.fq"; 
#my $file_un2=$inName.".un2.fq";  
#my $file_out=$inName.".stitched.fastq";
my $file_join=$inName.".assembled.fastq";
my $file_un1=$inName.".unassembled.forward.fastq"; 
my $file_un2=$inName.".unassembled.reverse.fastq";  
my $file_out=$inName.".stitched.fastq";

open (INjoin, "<$file_join") or die "Can't Open file: $file_join";
open (INun1,"<$file_un1") or die "Can't Open file: $file_un1";
open (INun2,"<$file_un2") or die "Can't Open file: $file_un2";
open (OUT,">$file_out") or die "Can't Open file: $file_out";

my $cc=0;
while (my $line = <INjoin>){
	$cc++; 
	if ($cc-4*int($cc/4)==1) {
            #print $line; <STDIN>; 
		chomp $line;
                $line =~ s/^@//;
                my @hdr = split(/\s+/,$line);
=pod
                my $info = $line;
                my @name1 = split(/\s+/,$line);  # get rid of the mate info 1:N:0:1
                my @name2 = split(/\:/, $name1[0]); 
                print "name2=",@name2,"\n";
                my $name = $name2[3];
                print "name=",$name,"\n";
		for(my $i=4; $i <= $#name2; $i++){ 
			$name=join("\:", $name, $name2[$i]); 
		}
                print OUT "\@join:$name\n";
 		print OUT "\@join:$line\n";
=cut
 		print OUT "\@join:$hdr[0]\n";
	}
	else { print OUT $line; }
}

$cc=0;
while (my $line = <INun1>){
	$cc++;
	if ($cc-4*int($cc/4)==1) {
		my $info = $line;
		chomp $line;
                $line =~ s/^@//;
                my @hdr = split(/\s+/,$line);
=pod
		my @name1 = split(/\s+/,$line);  # get rid of the mate info 1:N:0:1
		my @name2 = split(/\:/, $name1[0]); 
		my $name = $name2[3];
		for(my $i=4; $i <= $#name2; $i++){ 
			$name=join("\:", $name, $name2[$i]); 
		}
		print OUT "\@un1:$name\n";
		print OUT "\@un1:$line\n";
=cut
		print OUT "\@un1:$hdr[0]\n";
	}
	else { print OUT $line; }
}

$cc=0; 
while (my $line = <INun2>){
	$cc++;
	if ($cc-4*int($cc/4)==1) {
		chomp $line;
                $line =~ s/^@//;
                my @hdr = split(/\s+/,$line);
=pod
		my $info = $line;
		my @name1 = split(/\s+/,$line);  # get rid of the mate info 1:N:0:1
		my @name2 = split(/\:/, $name1[0]); 
		my $name = $name2[3];
		for(my $i=4; $i <= $#name2; $i++){ 
			$name=join("\:", $name, $name2[$i]); 
		}
		print OUT "\@un2:$name\n";
		print OUT "\@un2:$line\n";
=cut
		print OUT "\@un2:$hdr[0]\n";
	}
	else { print OUT $line; }
}

close INjoin; 
close INun2; 
close INun1; 

exit;
