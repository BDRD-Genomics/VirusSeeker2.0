# Steps in VirusSeeker

STAGE 1: Copy reads to output directory
STAGE 2: Quality control with fastp
STAGE 3a: Map reads to reference genome
STAGE 3b: Skip host removal, send output to QC
STAGE 3c: Map ONT reads to reference genome
STAGE 4a: Perform assembly using metaSPAdes: both paired end and single end reads in the same step
STAGE 4b: Perform long-read OR hybrid assembly using dragonflye
STAGE 4c: Perform hybrid assembly using unicycler
STAGE 4d: Perform assembly using SPAdes: both paired end and single end reads in the same step
STAGE 5a: Map reads back to contigs for short reads
STAGE 5b: Map reads back to contigs for long reads
STAGE 4b/5c: Skip assembly/contig mapping
STAGE 6: Stitch paired reads together if unmapped to assembly
STAGE 7: Combine sequences into single file
STAGE 8: This step will cluster sequences using mmseqs to reduce redundancy
STAGE 9: Split unique reads to small files for RepeatMasker
STAGE 10: Run RepeatMasker
STAGE 11: Check RepeatMasker output
STAGE 12: Split input file for blastn against viral-only database
STAGE 13: Run mmseqs2 against viral-only database
STAGE 14: Parse output of mmseqs against viral-only database
STAGE 15: Split input file for blastx against viral-only database
STAGE 16: Run blastx against viral-only database
STAGE 17: Parse output of blastx against viral-only database
STAGE 18: Split input file for megablast against NT
STAGE 19: Run megablast against NT
STAGE 20: Parse output for megablast against NT
STAGE 21: Split input file for blastx against NR
STAGE 22: Run diamond blastx against NR
STAGE 23: Parse output for blastx against NR
STAGE 24: Map reads to viral reference genome
STAGE 24B: Map reads to viral reference genome
STAGE 25: Generate Assignment Report
STAGE 26: Generate Assignment Summary
STAGE 27: Generate Phage Report
STAGE 28: Generate Phage Summary
STAGE 29: Generate Supplemental Outputs
