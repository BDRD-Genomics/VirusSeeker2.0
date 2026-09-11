---
title: 'VirusSeeker 2.0: A Network-Independent, Containerized Pipeline for the Discovery of Highly Divergent Viral Genomes'
tags:
  - bioinformatics
  - virus discovery
  - biosurveillance
  - de novo assembly
  - short-read-sequencing
  - long-read sequencing
authors:
- family-names: "Rice"
  given-names: "Gregory K"
  orcid: "https://orcid.org/0000-0001-8509-6379"
  affiliation: "1, 2"
- family-names: "Paskey"
  given-names: "Adrian C"
  orcid: "https://orcid.org/0000-0003-4575-3092"
- family-names: "Thomas"
  given-names: "Quinn K"
  orcid: "https://orcid.org/0009-0008-2920-4484"
  affiliation: "1, 2"
- family-names: "Long"
  given-names: "Kyle A."
  orcid: ""
  affiliation: "1, 2"
- family-names: "Cicalo"
  given-names: "Anthony"
  orcid: ""
  affiliation: "1, 2"
- family-names: "Philipson"
  given-names: "Casandra W."
  orcid: ""
  affiliation: "3"
- family-names: "Zhao"
  given-names: "Guoyan"
  orcid: "https://orcid.org/0000-0001-5615-6774"
  affiliation: "4, 5, 6, *"
- family-names: "Cer"
  given-names: "Regina Z."
  orcid: "https://orcid.org/0000-0002-2395-980X"
  affiliation: "1"
- family-names: "Bishop-Lilly"
  given-names: "Kimberly A."
  orcid: "https://orcid.org/0000-0002-5744-8493"
  affiliation: "1"
affiliations:
 - name: Genomics and Bioinformatics Department, Biological Defense Research Directorate, Naval Medical Research Command-Frederick, United States
   index: 1
 - name:  Leidos, USA
   index: 2
 - name:  Defense Threat Reduction Agency, USA
   index: 3
 - name:  Department of Genetics, Washington University School of Medicine, USA
   index: 4
 - name:  Department of Neurology, Washington University School of Medicine, USA
   index: 5
 - name:  Departments of Pathology and Immunology, Washington University School of Medicine, USA
   index: 6
 - name:  "*Deceased"
   index: "*"
date: 3 September 2026
bibliography: paper.bib
---

# Summary

`VirusSeeker 2.0` is a specialized bioinformatics pipeline engineered specifically for high-sensitivity detection and characterization of viral sequences within complex metagenomic datasets. Operating either as a standalone tool or as a complementary specialized tool to other broad-domain profiling tools, VirusSeeker 2.0 analyzes assembled contigs and stitched reads to detect low-abundance viral sequences within samples dominated by host or bacterial sequence. Similar to other pipelines, it automates quality control (BBDuk (1), fastp (2)), host sequence removal (BBMap (1), minimap2 (3)), and de novo assembly (metaSPAdes (4), SPAdes (5), Unicycler (6), Dragonflye (7)). However, it diverges in its subsequent taxonomic classification by mining not only assembled contigs but also unassembled viral reads that were not incorporated into the assembly using DIAMOND (8) and MEGAN (9). The primary output is a report cataloging viral reads and contigs detected through nucleotide and amino acid level sequence comparisons. This dual-level search strategy is optimized to uncover novel or highly divergent viral sequences that lack close reference matches. 

# Statement of need

`VirusSeeker 2.0` addresses a computational challenge in viral discovery from high throughput sequencing data: viral reads and contigs may be highly divergent from nearest sequenced neighbors and can therefore be missed by conventional reference-based analyses. Furthermore, comprehensive viral discovery requires multiple tools, many of which assume active internet connectivity and are therefore difficult to use in air-gapped high-performance computing (HPC) environments. VirusSeeker 2.0 overcomes these infrastructure limitations by packaging a reproducible viral discovery workflow in portable containers for offline execution.  

# State of the field                                                                                                                  

`VirusSeeker 2.0` provides high sensitivity viral detection and taxonomic classification from metagenomic data, purposely built for high performance in air gapped environments. To evaluate VirusSeeker 2.0’s capacity to detect challenging viral genomes using short-read and hybrid sequencing data, we benchmarked it against our complementary in-house pipeline, MetaDetector, and four widely used tools: geNomad v1.8.0 (30), Chan Zuckerberg ID (CZID v6.0 (31, 32), Kraken2 v2.1.3 (33), and Mash v2.3 (34).
Each tool selected for this benchmarking exercise represents a distinct methodological framework, and was executed using default parameters: 
  - **`VirusSeeker 2.0`** is a discovery pipeline combining host subtraction, assembly, and multi-tier BLAST classification with false positive filtering.
  - **geNomad** uses profile HMMs and marker-based methods to identify  viral and plasmid sequences.
  - **CZID** and **MetaDetector** perform de novo short-read assembly via SPAdes (16), followed by BLAST-based contig classification (35).
  - **Kraken2** uses an exact k-mer matching approach for classification of the lowest common ancestor (LCA).
  - **Mash** uses hash containment screening against RefSeq reference sketches.

Two datasets from NCBI’s Sequence Read Archive (SRA) were selected for benchmarking: 
1. SRR10179613, Illumina short reads generated from a fruit bat body swab, selected to test sensitivity in detecting a low-abundance, recently discovered virus (Dawn bat paramyxovirus, DbPV) against a background of host and bacterial reads.
2. PRJNA587334, matched Oxford Nanopore MinION and Illumina MiSeq read sets generated from an mpox virus (MPXV) clinical isolate (29), selected to assess hybrid long/short read assembly and classification.                                                                                         



## Benchmarking Results
[Table 1](#Table-1) summarizes the degree to which all evaluated tools identified the target species or a closely related near-neighbor. 
  - In the environmental bat-swab sample (Table 1A), VirusSeeker 2.0 showed the highest sensitivity, detecting 1,102 reads and assembling 36 contigs assigned to DbPV. In contrast, geNomad resolved DbPV only to the family level (Paramyxoviridae), while MetaDetector, CZID, and Kraken2 provided species-level resolution by identifying closely related near-neighbor lineages.
  - All evaluated tools successfully classified Mpox virus in the hybrid assembly dataset (Table 1B). VirusSeeker 2.0 assigned 11,814 reads to MPXV; however, its strict false-positive filtering excluded MPXV contigs that aligned to vaccinia virus synthetic constructs (e.g., AY965296.1). CZID and MetaDetector identified the highest number of MPXV reads, with 54,456 and 49,047 reads respectively. MetaDetector, CZID, and Kraken2 each produced comparable numbers of MPXV contigs (12, 13, and 9 contigs respectively).
  - These results demonstrate the  rigorous false-positive filtering of VirusSeeker 2.0 to optimize for the detection of low-abundance novel viruses amid high background, which is unique in comparison to the other tools that excel more broadly at assembly and contig recovery for all domains of life and therefore these results exemplify the specific use case for which VirusSeeker 2.0 was developed – virus discovery.

<a id="Table-1"></a>
<i> **Table 1.** Summary of benchmarking results via VirusSeeker 2.0, MetaDetector, geNomad, CZID, kraken2, and mash.</i> \
1A. Detection of virus
| Tool | Specificity of DbPV detection | Specific Assignment | # reads classified at the lowest assignment | # contigs classified at the lowest assignment |
| --- | --- | --- | --- | --- |
| VirusSeeker 2.0 | high | DbPV | 1102 | 36 |
| MetaDetector | moderate | Mumps/bat mumps orthorubulavirus (TaxIDs: 2560602; 2560195) | 63 | 1 |
| geNomad | low | Paramyxoviridae (TaxID: 11158) | n/a | n/a |
| CZID | moderate | Bat mumps orthorubulavirus (TaxID: 2560340) | 16 | 0 |
| kraken2 | moderate | Bat Paramyxovirus Epo_spe/AR1/DRC/2009 (NC_038271) | 2 | 1 |
| mash | moderate | Bat Paramyxovirus Epo_spe/AR1/DRC/2009n(NC_038271) | n/a | n/a |

1B. Handling of hybrid mpox virus data
| Tool | Specificity of DbPV detection | Specific Assignment | # reads classified at the lowest assignment | # contigs classified at the lowest assignment |
| --- | --- | --- | --- | --- |
| VirusSeeker 2.0* | high | MPXV (TaxID: 10244) | 11814 | 0* |
| MetaDetector | high | MPXV (TaxID: 10244) | 49047 | 12 | 
| geNomad | low | Poxviridae (TaxID: 10240) | n/a | n/a |
| CZID** | high | MPXV (TaxID: 10244) | 54456 | 13 |
| kraken2 | high | MPXV (TaxID: 10244) | 35595 | 9 |
| mash | high | MPXV (NC_063383) | n/a | n/a |

*VirusSeeker 2.0 contigs aligned to synthetic construct clone records (e.g., AY965296.1) and were filtered out by the pipeline's non-viral false-positive filter, demonstrating the benefit of running VirusSeeker 2.0 alongside MetaDetector.\
**CZ ID was run locally; results reflect short-read input only, as hybrid input is not supported.
 

# Software design

`VirusSeeker 2.0` 0 modernizes the original BLAST-based VirusSeeker framework while retaining its emphasis on detecting low-prevalence, highly divergent eukaryotic viruses while enforcing rigorous false-positive filtration. The redesigned workflow addresses throughput, memory, and database bottlenecks of legacy metagenomic discovery pipelines via the following elements: 
  - Dual analysis tracks (previously Discovery vs. Virome) were unified by replacing disparate Perl scripts with a single-entry point, therefore streamlining execution using Slurm and containerized (Docker) HPC environments. 
  - Post-assembly PEAR read stitching and back-mapping to contigs replaced the legacy pipeline pre-assembly stitching to maximize detection sensitivity of unassembled or divergent viral fragments.
  - Clustering and repeat masking were implemented via MMseqs2 to mitigate the bottleneck of RepeatMasker and BLAST. Sequence sets are partitioned across available threads to enable high throughput screening (500 chunks, >36k sequences).
  - Multi-tiered alignment using MMseqs2 nucleotide search and DIAMOND BLASTX against virus-specific and full NT/NR databases improves compute time over traditional BLASTX searching against NCBI NR by taking advantage of DIAMOND’s algorithm without compromising sensitivity for novel viral homologs.
  - The development of an automated Accession.Version tracking and normalized length-adjusted read counts provides abundance quantification that corrects for viral genome size disparities by generating reads per million (RPM) and family-level normalized abundance reports alongside strict non-viral filtering.

# Research impact statement

Methods for detecting novel and/or low-abundance pathogens amidst background host genomes typically overwhelm compute infrastructure due to the best practice of querying distant protein homologs. By modernizing legacy GI-based parsing to Accession.Version schemas, moving read stitching to post-assembly, and updating to accelerated alignment heuristics (MMseqs2/DIAMOND), VirusSeeker 2.0 has enabled sensitive virus discovery in an air gapped environment for five peer reviewed publications (10-14). VirusSeeker 2.0 was implemented using Docker (45) to provide reach-back support and train 19 personnel at a U.S. Department of Defense laboratory overseas. Trainees were able to individually analyze a sample from long and short read data through identifying the organism(s) of interest and performed advanced characterization, including phylogenetic analysis, of the target organism(s).

# AI usage disclosure

Generative AI was not used in any part of the software creation or documentation, which was designed and prepared by human authors. GenAI.mil was used to proofread sections of the written documentation and any edits resulting from those recommendations were vetted and implemented by human authors.

# Acknowledgements
This work was supported by Navy WUN A1417 and Global Emerging Infections Surveillance (GEIS) Branch ProMIS ID P0054_23_NM to KAB-L.
The views expressed in this article are those of the authors and do not necessarily reflect the official policy or position of the Department of the Navy, Department of Defense, nor the U.S. Government. Some authors are employees of the U.S. Government. This work was prepared as part of their official duties. Title 17 U.S.C. §105 provides that “Copyright protection under this title is not available for any work of the United States Government”. Title 17 U.S.C. §101 defines a U.S. Government work as a work prepared by a military service member or employee of the U.S. Government as part of that person’s official duties.

We dedicate this work to the memory of Dr. Guoyan Zhao, Associate Professor of Genetics and Neurology at Washington University School of Medicine.  As the primary creator and lead developer of the original VirusSeeker pipeline, her pioneering work laid the critical computational foundation upon which `VirusSeeker 2.0` was built.

# References
1.	Bushnell B. 2014. BBMap: a fast, accurate, splice-aware aligner. Lawrence Berkeley National Lab (LBNL).
2.	Chen S. 2023. Ultrafast one‐pass FASTQ data preprocessing, quality control, and deduplication using fastp. Imeta 2:e107.
3.	Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100.
4.	Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. 2017. metaSPAdes: a new versatile metagenomic assembler. Genome Res 27:824-834.
5.	Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, Lesin VM, Nikolenko SI, Pham S, Prjibelski AD. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology 19:455-477.
6.	Wick RR, Judd LM, Gorrie CL, Holt KE. 2017. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS computational biology 13:e1005595.
7.	Petit III RA. 2024. Dragonflye: assemble bacterial isolate genomes from nanopore reads. Github2021.
8.	Buchfink B, Reuter K, Drost HG. 2021. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat Methods 18:366-368.
9.	Huson DH, Auch AF, Qi J, Schuster SC. 2007. MEGAN analysis of metagenomic data. Genome research 17:377-386.
10.	Adhikari BN, Paskey AC, Frey KG, Bennett AJ, Long KA, Kuhn JH, Hamilton T, Glang L, Cer RZ, Goldberg TL, Bishop-Lilly KA. 2024. Virome profiling of fig wasps (Ceratosolen spp.) reveals virus diversity spanning four realms. Virology 591:109992.
11.	Bennett AJ, Paskey AC, Kuhn JH, Bishop-Lilly KA, Goldberg TL. 2020. Diversity, Transmission, and Cophylogeny of Ledanteviruses (Rhabdoviridae: Ledantevirus) and Nycteribiid Bat Flies Parasitizing Angolan Soft-Furred Fruit Bats in Bundibugyo District, Uganda. Microorganisms 8.
12.	Paskey AC, Lim XF, Ng JHJ, Rice GK, Chia WN, Philipson CW, Foo R, Cer RZ, Long KA, Lueder MR, Glang L, Frey KG, Hamilton T, Mendenhall IH, Smith GJ, Anderson DE, Wang LF, Bishop-Lilly KA. 2023. Genomic Characterization of a Relative of Mumps Virus in Lesser Dawn Bats of Southeast Asia. Viruses 15.
13.	Paskey AC, Ng JHJ, Rice GK, Chia WN, Philipson CW, Foo RJH, Cer RZ, Long KA, Lueder MR, Frey KG, Hamilton T, Mendenhall IH, Smith GJ, Wang LF, Bishop-Lilly KA. 2020. The temporal RNA virome patterns of a lesser dawn bat (Eonycteris spelaea) colony revealed by deep sequencing. Virus Evol 6:veaa017.
14.	Bennett AJ, Paskey AC, Ebinger A, Pfaff F, Priemer G, Hoper D, Breithaupt A, Heuser E, Ulrich RG, Kuhn JH, Bishop-Lilly KA, Beer M, Goldberg TL. 2020. Relatives of rubella virus in diverse mammals. Nature 586:424-428.
