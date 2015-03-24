# preprocess_app.py
import sys

"""
Preprocessing of fastq sequence files for a single sample,
in either paired format or single end format and potentially both

Preprocessing involves the following applications
1. Contaminant screening: using Bowtie local very-sensitive
2. de-duplications: Super_Deduper
3. Quality trimming, polyA/T trimming: Sickle2
    File handle in -> STDOUT
4. Overlapping paired-end reads, Flash2
5. Normalization, kmer_filter from Stacks

TODO: 
Modify each of the applications above to accept a tabbed sequence format
to facilitate fast streaming of reads STDIN, STDOUT and moduleration

tabbed sequence format is the same as described by FLASH
" In this mode you should provide a single input file,
 each line of which must contain either a read pair (5 fields)
 or a single read (3 fields)."

This format will facilitate processing both single and paired-end read seamlessly.
"""


class preprocessCMD:
    """
    preprocessing application for expHTS experiments
    """
    def __init__(self):
        pass

    def execute(self, args):
        # ----------------------- options input files -----------------------
        if args.samples_file is None:
            sFile = None
            sys.stderr.write("No sample file identified\n")
        else:
            sFile = args.samples_file
        print sFile

        return 0
