#!/usr/bin/env bash
set -Eeuo pipefail

SRA_ACCESSION="SRR22470068"
SARS_NT="NC_045512.2"
HUMAN_MT="NC_012920.1"
ALU_SX="U14574.1"

src=/database/test_source
inp=/input
mkdir -p "$src" "$inp" /scratch/build/sra

EFETCH='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

fetch_fasta() {
    local db="$1" id="$2" out="$3"
    if [[ -s "$out" ]]; then
        echo "Reusing $out"
        return 0
    fi
    curl -fL --retry 5 --retry-delay 2 \
      "$EFETCH?db=$db&id=$id&rettype=fasta&retmode=text" \
      -o "$out"
    grep -q '^>' "$out"
}

# Real SARS-CoV-2 reference genome.
fetch_fasta nuccore "$SARS_NT" "$src/NC_045512.2.fa"

fetch_fasta nuccore "$HUMAN_MT" "$src/NC_012920.1.fa"

# Real human Alu-Sx consensus sequence from GenBank for the small RepeatMasker
# custom library used by the functional test.
fetch_fasta nuccore "$ALU_SX" "$src/U14574.1_AluSx.fa"

# Real RefSeq proteins encoded by NC_045512.2.
protein_ids='YP_009724389.1,YP_009724390.1,YP_009724391.1,YP_009724392.1,YP_009724393.1,YP_009724394.1,YP_009724395.1,YP_009724396.1,YP_009724397.2,YP_009725255.1,YP_009725295.1'
if [[ ! -s "$src/NC_045512.2.proteins.fa" ]]; then
    curl -fL --retry 5 --retry-delay 2 \
      "$EFETCH?db=protein&id=$protein_ids&rettype=fasta&retmode=text" \
      -o "$src/NC_045512.2.proteins.fa"
fi
grep -q '^>' "$src/NC_045512.2.proteins.fa"


R1_OUT="$inp/${SRA_ACCESSION}_R1.fastq.gz"
R2_OUT="$inp/${SRA_ACCESSION}_R2.fastq.gz"

if [[ ! -s "$R1_OUT" || ! -s "$R2_OUT" ]]; then
    rm -f "$R1_OUT" "$R2_OUT" \
          "$inp/${SRA_ACCESSION}_1.fastq.gz" \
          "$inp/${SRA_ACCESSION}_2.fastq.gz"
    rm -rf "/scratch/build/sra/$SRA_ACCESSION"
    mkdir -p "/scratch/build/sra/$SRA_ACCESSION"

    cd "/scratch/build/sra/$SRA_ACCESSION"
    /opt/conda/envs/vs/bin/fasterq-dump \
      --split-files \
      --threads "${VS_TEST_DOWNLOAD_THREADS:-4}" \
      --outdir . \
      "$SRA_ACCESSION"

    test -s "${SRA_ACCESSION}_1.fastq"
    test -s "${SRA_ACCESSION}_2.fastq"

    pigz -p "${VS_TEST_DOWNLOAD_THREADS:-4}" -c \
      "${SRA_ACCESSION}_1.fastq" > "$R1_OUT"
    pigz -p "${VS_TEST_DOWNLOAD_THREADS:-4}" -c \
      "${SRA_ACCESSION}_2.fastq" > "$R2_OUT"

    test -s "$R1_OUT"
    test -s "$R2_OUT"
fi

cat > "$inp/README.test-data.txt" <<EOF_README
VirusSeeker public toy/test data set

All biological sequence data in the test profile are public biological records; no synthetic reads or sequences are generated.

Public read set:
  NCBI SRA: $SRA_ACCESSION
  Organism: Severe acute respiratory syndrome coronavirus 2
  Layout: paired-end Illumina

Public references used to build the compact test databases:
  SARS-CoV-2 genome: NC_045512.2
  SARS-CoV-2 RefSeq proteins: YP_009724389.1 ... YP_009725295.1
  Minimal human host reference: NC_012920.1 (human mitochondrial genome)
  RepeatMasker test library: U14574.1 (human Alu-Sx consensus)

The test profile is an execution/reproducibility test, not a sensitivity benchmark.
EOF_README

printf 'Downloaded real public test assets successfully. SRA=%s\n' "$SRA_ACCESSION"
