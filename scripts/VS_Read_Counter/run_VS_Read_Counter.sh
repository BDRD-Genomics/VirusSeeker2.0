###############################################################################
# Author: BDRD <usn.detrick.nmrc.mbx.genomics-reach-back@health.mil>
#
# License:
# VirusSeeker 2.0 is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3 of the License or any later version. For 
# details, please refer to https://www.gnu.org/licenses/
###############################################################################

#echo "Creating directories..."
mkdir -p $1/Supplemental_Outputs
mkdir -p $1/Supplemental_Outputs/Fungi
mkdir -p $1/Supplemental_Outputs/unassigned
mkdir -p $1/Supplemental_Outputs/Bacteria
mkdir -p $1/Supplemental_Outputs/Homo
mkdir -p $1/Supplemental_Outputs/Mus
mkdir -p $1/Supplemental_Outputs/other
mkdir -p $1/Supplemental_Outputs/Phage
mkdir -p $1/Supplemental_Outputs/Virus
mkdir -p $1/Supplemental_Outputs/Ambiguous

#echo "Generating extra assignment reports and read lists..."
for reportScript in $2/VS_Read_Counter/*.pl
do
	perl $reportScript $1;
done




#echo "Moving files..."
mv $1/*FungiReadsAssignmentReport $1/Supplemental_Outputs/Fungi
mv $1/*unassignedReadsAssignmentReport $1/Supplemental_Outputs/unassigned
mv $1/*BacteriaReadsAssignmentReport $1/Supplemental_Outputs/Bacteria
mv $1/*HomoReadsAssignmentReport $1/Supplemental_Outputs/Homo
mv $1/*MusReadsAssignmentReport $1/Supplemental_Outputs/Mus
mv $1/*otherReadsAssignmentReport $1/Supplemental_Outputs/other
mv $1/*PhageReadsAssignmentReport $1/Supplemental_Outputs/Phage
cp $1/*.AssignmentReport $1/Supplemental_Outputs/Virus
cp $1/*AmbiguousReadsAssignmentReport $1/Supplemental_Outputs/Ambiguous


mv $1/*FungiReads_all.fa $1/Supplemental_Outputs/Fungi
mv $1/*unassignedReads_all.fa $1/Supplemental_Outputs/unassigned
mv $1/*BacteriaReads_all.fa $1/Supplemental_Outputs/Bacteria
mv $1/*HomoReads_all.fa $1/Supplemental_Outputs/Homo
mv $1/*MusReads_all.fa $1/Supplemental_Outputs/Mus
mv $1/*otherReads_all.fa $1/Supplemental_Outputs/other
cp $1/*PhageReads_All.fa $1/Supplemental_Outputs/Phage
cp $1/*ViralReads_all.fa $1/Supplemental_Outputs/Virus
cp $1/*AmbiguousReads_all.fa  $1/Supplemental_Outputs/Ambiguous



#echo "Calculating accurate read counts for all 9 VirusSeeker classifications..."


for organism in $1/Supplemental_Outputs/*;
do
	sample=$(basename $1)
	python $2/VS_Read_Counter/VS_Read_Counter.py $1 $organism/*AssignmentReport $organism/${sample}_AccurateReadCounts.tsv
        organism_bn=$(basename $organism)
        if [[ $organism_bn = 'Virus' ]]
        then
        #    echo "Testing 123!"
                python $2/VS_Read_Counter/normalize_reads.py -i $organism/${sample}_AccurateReadCounts.tsv
        fi;
done
#echo "Done!"
