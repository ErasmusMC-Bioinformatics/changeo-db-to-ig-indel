# changeo-db-to-ig-indel
Script that creates valid Ig-Indel input from a Change-O DB and fastq file and runs Ig-Indel on it.

    usage: changeo-db-to-ig-indel.py [-h] -cdb CHANGE_O_DB --fastq FASTQ -o
                                     OUTPUT_DIR [--mafft MAFFT]
                                     [--ig-indel_jar IG_INDEL_JAR]
                                     [--ig-indel-num-seqs IG_INDEL_NUM_SEQS]
                                     [--ig-indel-hom-length IG_INDEL_HOM_LENGTH]
                                     [--ig-indel-min-quality IG_INDEL_MIN_QUALITY]

Dependencies
--------------------
[Python 2.7](https://www.python.org/)  
[Biopython](http://biopython.org/wiki/Biopython)  
