# Downloading Databases

## Set Paths
```
vsdir=/path/to/databases  
scripts=/path/to/VirusSeeker_2.0/scripts
```

## Download NR and NT Commands
```
update_blastdb.pl --decompress nr && mv nr $vsdir/nr_$(date +%d%b%y) 
ln -sf $vsdir/nr_$(date +%d%b%y) $vsdir/nr  
update_blastdb.pl --decompress core_nt && mv core_nt $vsdir/core_nt_$(date +%d%b%y)  
ln -sf $vsdir/core_nt_$(date +%d%b%y) $vsdir/core_nt 
```

## Download human reference genome
```
cd $vsdir/ref_genomes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

## Download taxdump and use taxonkit to extract Virus taxidlist
```
wget -cv -O $vsdir/taxdump_$(date +%d%b%y) https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf $vsdir/taxdump_$(date +%d%b%y)
ln -sf $vsdir/taxdump_$(date +%d%b%y) $vsdir/taxdump
taxonkit list --ids 10239 -I "" >  $vsdir/virus.taxid.txt
 ```

## Clustering VirusDBNR and VirusDBNT Database
```
mkdir $vsdir/VirusDBNR_$(date +%d%b%y) && ln -sf $vsdir/VirusDBNR_$(date +%d%b%y) $vsdir/VirusDBNR && cd $vsdir/VirusDBNR
blastdbcmd -taxidlist $vsdir/virus.taxid.txt -dbtype prot -db $vsdir/nr/nr > virus_nr.fasta 
mmseqs easy-linclust $vsdir/VirusDBNR/virus_nr.fasta $vsdir/VirusDBNR/virus_nr.clustr98_98 /tmp --min-seq-id 0.98 -c 0.98 --threads 255
```
```
mkdir $vsdir/VirusDBNT_$(date +%d%b%y) && ln -sf $vsdir/VirusDBNT_$(date +%d%b%y) $vsdir/VirusDBNT && cd $vsdir/VirusDBNT
blastdbcmd -taxidlist $vsdir/virus.taxid.txt -dbtype nucl -db $vsdir/nt/nt > virus_nt.fasta
mmseqs easy-linclust $vsdir/VirusDBNT/virus_nt.fasta $vsdir/VirusDBNT/virus_nt.clustr98_98 /tmp --split-memory-limit 500G --min-seq-id 0.98 -c 0.98 --threads 255
```

## Creating BLAST Database for VirusDBNT
```
cd $vsdir/VirusDBNT
makeblastdb -in virus_nt.clustr98_98_rep_seq.fasta -blastdb_version 5 -dbtype nucl -parse_seqids -out virus_nt
dustmasker -in virus_nt -infmt blastdb -parse_seqids -outfmt maskinfo_asn1_bin -out virus_nt.asnb
makeblastdb -in virus_nt -input_type blastdb -blastdb_version 5 -dbtype nucl -parse_seqids -mask_data virus_nt.asnb -out virus_nt -title "Nucleotide database built from viral taxonomy ids extracted from NCBI nt clustered with mmseqs2 and masked with dustmasker"
```

## Creating Diamond Databases for VirusDBNR and NR 
```
cd $vsdir/VirusDBNR
diamond makedb --in virus_nr.clustr98_98_rep_seq.fasta --db virus_nr.dmnd
mkdir $vsdir/nr_dmnd_$(date +%d%b%y) && ln -sf $vsdir/nr_dmnd_$(date +%d%b%y) $vsdir/nr_dmnd && cd $vsdir/nr_dmnd
blastdbcmd blastdbcmd -db $vsdir/nr/nr -entry all | pigz -9 > $vsdir/nr_dmnd/nr.gz
cd $vsdir/nr_dmnd
diamond makedb --in nr.gz --db nr.dmnd
```

## Create vhunter_acc.db
```
wget -cv -O $vsdir/taxdump/nucl_gb.accession2taxid.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget -cv -O $vsdir/taxdump/prot.accession2taxid.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
pigz -dc nucl_gb.accession2taxid.gz > nucl_gb.accession2taxid
pigz -dc nucl_gb.accession2taxid.gz > prot.accession2taxid
cd $vsdir/taxdump/
bash $scripts/build_db_acc_taxid_nucl.sh  
bash $scripts/build_db_acc_taxid_prot.sh
```