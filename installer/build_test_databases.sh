#!/usr/bin/env bash
set -Eeuo pipefail

src=/database/test_source

rm -rf \
  /database/ref /database/nr /database/nt /database/core_nt \
  /database/nr_dmnd /database/core_nt_mmseqs \
  /database/VirusDBNR /database/VirusDBNT /database/test

mkdir -p \
  /database/ref /database/nr /database/nt /database/core_nt \
  /database/nr_dmnd /database/core_nt_mmseqs \
  /database/VirusDBNR /database/VirusDBNT /database/test \
  /scratch/mmseqs/test_core /scratch/mmseqs/test_virus

# Real NCBI Alu-Sx consensus used as a compact RepeatMasker library.
cp "$src/U14574.1_AluSx.fa" /database/test/repeat_library.fa

# Minimal real human host reference.  Keep the historical runtime prefix expected
# by runVS.pl while documenting that test mode uses only the mitochondrial genome.
cp "$src/NC_012920.1.fa" /database/ref/GRCh38.p14
bowtie2-build /database/ref/GRCh38.p14 /database/ref/GRCh38.p14

# Real SARS-CoV-2 representative genome used by the viral-reference mapping step.
cp "$src/NC_045512.2.fa" /database/ref/ref_viruses_rep_genomes
bowtie2-build /database/ref/ref_viruses_rep_genomes /database/ref/ref_viruses_rep_genomes

# Tiny but fully real protein databases.
makeblastdb \
  -in "$src/NC_045512.2.proteins.fa" \
  -dbtype prot -parse_seqids -out /database/nr/nr

diamond makedb \
  --in "$src/NC_045512.2.proteins.fa" \
  --db /database/nr_dmnd/nr

diamond makedb \
  --in "$src/NC_045512.2.proteins.fa" \
  --db /database/VirusDBNR/virus_nr

# Tiny real nucleotide background: SARS-CoV-2 plus a small human sequence.
cat "$src/NC_045512.2.fa" "$src/NC_012920.1.fa" > "$src/core_nt_test.fa"
makeblastdb \
  -in "$src/core_nt_test.fa" \
  -dbtype nucl -parse_seqids -out /database/nt/nt
makeblastdb \
  -in "$src/core_nt_test.fa" \
  -dbtype nucl -parse_seqids -out /database/core_nt/core_nt

mmseqs createdb "$src/core_nt_test.fa" /database/core_nt_mmseqs/core_nt
mmseqs createdb "$src/NC_045512.2.fa" /database/VirusDBNT/virus_nt


diamond blastx \
  --query "$src/NC_045512.2.fa" \
  --db /database/VirusDBNR/virus_nr.dmnd \
  --out "$src/test_blastx_preflight.tsv" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 5 --threads 1

test -s "$src/test_blastx_preflight.tsv"
head -n 3 "$src/test_blastx_preflight.tsv" >&2

test -s /database/nr_dmnd/nr.dmnd
test -s /database/VirusDBNR/virus_nr.dmnd
test -s /database/core_nt_mmseqs/core_nt.dbtype
mmseqs dbtype /database/core_nt_mmseqs/core_nt >/dev/null
test -s /database/VirusDBNT/virus_nt.dbtype
mmseqs dbtype /database/VirusDBNT/virus_nt >/dev/null
test -s /input/SRR22470068_R1.fastq.gz
test -s /input/SRR22470068_R2.fastq.gz
