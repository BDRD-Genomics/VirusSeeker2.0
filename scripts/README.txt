In this version:

Modified BLASTn and BLASTx parser to bin sequences that hit both
virus and other species to "Ambiguous" bin. Sequences only hit virus with
significant e value are counted as true "viral".
Finished modification 2014/02/04.
Updated to keep sequences that hit both virus and vector sequences.

Use Newbler for 1st step assembly and 2nd step singleton assembly
Use Phrap for 2nd step contig assembly.

Added BLASTx virus only database step.

Modified BLASTn nt, BLASTx NR and BLASTX VIRUSDB parser to generate
best hit alignment information in the parsed file. 

BLASTn VIRUSDB, BLASTX VIRUSDB, MegaBLAST nt, BLASTn nt, BLASTx NR

Used CD-HIT 95% cutoff to cluster reads and use respresentative + 
top 3 longest reads for assembly as the worm data suggest this is 
a better strategy.
