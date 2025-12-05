#!/bin/bash
indir=$1
outfile=$2
find ${indir} \( -name "*.AssignmentSummary" -o -wholename "*Supplemental_Outputs/Virus/*" \) | tar -cT - | pigz -9 -p 127 > ${outfile}
