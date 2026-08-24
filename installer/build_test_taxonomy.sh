#!/usr/bin/env bash
set -Eeuo pipefail

release_dir="${1:?usage: build_test_taxonomy.sh RELEASE_DIR}"
cd "$release_dir"

rm -f taxdump.tar.gz.part
if [[ -s taxdump.tar.gz ]] && gzip -t taxdump.tar.gz 2>/dev/null; then
    echo 'Existing taxdump.tar.gz passed gzip validation; reusing it.'
else
    rm -f taxdump.tar.gz
    wget --tries=5 --timeout=30 \
      -O taxdump.tar.gz.part \
      https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    gzip -t taxdump.tar.gz.part
    mv -f taxdump.tar.gz.part taxdump.tar.gz
fi

gzip -t taxdump.tar.gz
tar -xzf taxdump.tar.gz

printf 'accession\taccession.version\ttaxid\tgi\n' > nucl_gb.accession2taxid
printf 'NC_045512\tNC_045512.2\t2697049\t0\n' >> nucl_gb.accession2taxid
printf 'NC_012920\tNC_012920.1\t9606\t0\n' >> nucl_gb.accession2taxid

printf 'accession\taccession.version\ttaxid\tgi\n' > prot.accession2taxid
for acc in \
  YP_009724389.1 YP_009724390.1 YP_009724391.1 YP_009724392.1 \
  YP_009724393.1 YP_009724394.1 YP_009724395.1 YP_009724396.1 \
  YP_009724397.2 YP_009725255.1 YP_009725295.1
do
    printf '%s\t%s\t2697049\t0\n' "${acc%.*}" "$acc" >> prot.accession2taxid
done

/opt/conda/envs/vs/bin/python /installer/build_test_virus_taxids.py \
  nodes.dmp virus.taxid.txt 10239

test -s virus.taxid.txt
