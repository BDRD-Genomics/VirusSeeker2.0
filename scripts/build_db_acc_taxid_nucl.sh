###############################################################################
# Author: Kyle Long <kyle.a.long8.ctr@mail.mil>
# Additional Contact: Regina Cer <regina.z.cer.civ@mail.mil>
#
# License:
# VirusSeeker 2.0 is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, version 3 of the License or any later version. For 
# details, please refer to https://www.gnu.org/licenses/
###############################################################################


#!/bin/bash -x

echo "CREATE TABLE acc_taxid_nucl (accession text, accession_version text PRIMARY KEY, tax_id integer, gi integer);" | sqlite3 vhunter_acc.db
echo -e ".separator \"\t\"\n.import ./nucl_gb.accession2taxid acc_taxid_nucl\n"| sqlite3 vhunter_acc.db
