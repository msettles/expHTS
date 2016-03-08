# preprocess_app.py
import argparse
import signal
import os
from preprocessCMD import preprocessCMD
from mappingCMD import mappingCMD
from htseqcountCMD import htseqCMD
from spadesCMD import spadesCMD
from forcepairCMD import forcepairCMD
from kmerFilterCMD import kmerFilterCMD
version_num = "1.0"
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

# get process group number to know what to kill
# get process group number to know what to kill
pgid = os.getpgid(os.getpid())


def signal_handler(signal, frame):

    # the signal handler to kill the entire process group
    # incase Cntr + C is sent (sub processes aren't killed
    print "Cntr + c was hit - ending process group number " + str(pgid)
    import glob
    for f in glob.glob(".screening_cont*"):
        os.remove(f)
    os.killpg(pgid, 9)


def htseqParser(subparser):
    htseq_parser = subparser.add_parser('htseq', help='runs htseq-count on mapped (bam) files')
    htseq_parser.add_argument('-f', '--samplesfile', help='The filename of the sample file [default: %(default)s]', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt')
    htseq_parser.add_argument('-r', '--bamFolder', help='Directory where the sequence data is stored [default: %(default)s]', action='store', type=str, dest='readFolder', metavar='FOLDER', default='03-BWA')
    htseq_parser.add_argument('-R', "--referenceGTF", help='Reference gtf to count against', action='store', type=str, dest='refGTF', metavar='REFERENCE GFT', default='')
    # htseq_parser.add_argument('-o', "--order", help='pos or name - [default name]', action='store', type=str, dest='order', metavar='ORDER', default='name');
    htseq_parser.add_argument('-s', "--stranded", help='yes, no, or reverse - [default: %(default)s]', action='store', type=str, dest='stranded', metavar='STRANDED', default='yes')
    htseq_parser.add_argument('-a', "--minaqual", help='skip all reads with alignment quality lower than the given minimum value [default: %(default)s]', action='store', type=str, dest='minQual', metavar='MINAQUAL', default='10')
    htseq_parser.add_argument('-y', "--type", help='feature type (3rd column in GFF file) to be used, all features of other type are ignored [default suitable for ensembly GTF files: %(default)s]', action='store', type=str, dest='type', metavar='TYPE', default='exon')
    htseq_parser.add_argument('-i', "--idattr", help='GFF attribute to be used as feature ID [default: %(default)s]"', action='store', type=str, dest='idattr', metavar='IDATTR', default='gene_id')
    htseq_parser.add_argument('-m', "--mode", help='union, intersection-strict, intersection-nonempty - [default: %(default)s]', action='store', type=str, dest='mode', metavar='MODE', default='union')
    htseq_parser.add_argument('-F', '--final-folder', help='folder name in which the sequences will go [default: %(default)s]', action='store', type=str, default="04-HTseqCounts", dest="finalDir", metavar='DIRECTORY')
    htseq_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default: %(default)s]', action='store_true', dest='force', default=False)

    return htseq_parser


def mappingParser(subparser):
    mapping_parser = subparser.add_parser('mapping', help='maps reads to a reference sequence and post processes')
    mapping_parser.add_argument('-f', '--samplesfile', help='The filename of the sample file [default: %(default)s]', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt')
    mapping_parser.add_argument('-r', '--readFolder', help='Directory where the sequence data is stored [default: %(default)s]', action='store', type=str, dest='readFolder', metavar='FOLDER', default='02-Cleaned')
    mapping_parser.add_argument('-R', "--reference", help='Reference fasta to map against', action='store', type=str, dest='refFasta', metavar='REFERENCE', default='')
    mapping_parser.add_argument('-i', '--force-index', help='overwrites old index files [default: %(default)s]', action='store_true', dest='forceIndex', default=False)
    mapping_parser.add_argument('-M', "--mappingAlgorithm", help='Mapping algorithm bwa or bowtie2 [default: %(default)s]', action='store', type=str, dest='mapping', metavar='ALGORITHM', default='bwa')
    mapping_parser.add_argument('-n', "--sortByReadID", help="When sorting bam files, sort by read ID (samtools -n option), for compatability with htseq-count [default: %(default)s]", action='store_true', dest='sortByReadID', default=False)
    mapping_parser.add_argument('-s', "--ignoreSingles", help="Ignore any single-end files, for compatability with htseq-count [default: %(default)s]", action='store_true', dest="ignoreSingles", default=False)
    mapping_parser.add_argument('-F', '--final-folder', help='folder name in which the sequences will go [default: %(default)s]', action='store', type=str, default="03-BWA", dest="finalDir", metavar='DIRECTORY')
    mapping_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default: %(default)s]', action='store_true', dest='force', default=False)
    mapping_parser.add_argument('-S', '--forcePairs', help='Force pairs [default: %(default)s]', action='store_true', dest='forcePairs', default=False)
    mapping_parser.add_argument('-t', '--threads', help='Number of threads to be used [default: %(default)s]', action='store', type=str, dest='threads', metavar='THREADS', default='20')

    return mapping_parser


def preprocessParser(subparser):
    expHTS_parser = subparser.add_parser('preprocess', help='runs the expHTS preprocessing pipeline')
    expHTS_parser.add_argument('-f', '--samplesfile', help='The filename of the sample file [default: %(default)s]', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt')
    expHTS_parser.add_argument('-S', '--forceSplit', help='Forces splits of SE reads [default: %(default)s]', action='store_true', dest='forceSplit', default=False)
    expHTS_parser.add_argument('-A', '--adapterfasta', help='folder name with adapter sequences in fasta format [default: %(default)s]', action='store', type=str, default=r'<(printf ">TruSeq_forward_contam\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[NNNNNN]ATCTCGTATGCCGTCTTCTGCTTGAAAAA\n>TruSeq_reverse_contam\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA")',  dest='adapter', metavar='CONTAMNANTS-FOLDER')
    expHTS_parser.add_argument('-d', '--directory', help='Directory where the raw sequence data is stored [default: %(default)s]', action='store', type=str, dest='samplesDirectory', metavar='DIRECTORY', default='00-RawData')
    expHTS_parser.add_argument('-q', '--quality', help='Quality score to use during lucy trimming [default: %(default)s]', action='store', type=str, dest='qualityTrim', metavar='QUALITY', default='20')
    expHTS_parser.add_argument('-m', '--miniumumLength', help='Discard reads less than minimum length [default: %(default)s]', action='store', type=str, dest='minLength', metavar='MINIUMUMLENGTH', default='50')
    expHTS_parser.add_argument('-o', '--overlap', help='Overlap parameter for flash [default: %(default)s]', action='store', type=str, dest='overlapFlash', metavar='OVERLAP', default='700')
    expHTS_parser.add_argument('-O', '--skip-overlap', help='do not perform the overlapping using flash [default: %(default)s]', action='store_true', dest='skipFlash',  default=False)
    expHTS_parser.add_argument('-s', '--skip-duplicates', help='do not preform the deduplication step [default: %(default)s]', action='store_true',  dest='skipDup',  default=False)
    expHTS_parser.add_argument('-c', '--contaminates-folder', help='folder name with contaminate sequences in fasta format [default: %(default)s]', action='store', type=str, default='',  dest='contaminantsFolder', metavar='CONTAMNANTS-FOLDER')
    expHTS_parser.add_argument('-a', '--polyA', help='perform polyA trimming in sickle [default: %(default)s]', action='store_true', dest='polyTrim', default=False)
    expHTS_parser.add_argument('-F', '--final-folder', help='folder name in which the sequences will go [default: %(default)s', action='store', type=str, default="02-Cleaned", dest="finalDir", metavar='DIRECTORY')
    expHTS_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default: %(default)s]', action='store_true', dest='force', default=False)
    expHTS_parser.add_argument('-t', '--threads', help='Threads for bowtie2 [default: %(default)s]', action='store', type=str, dest='threads', metavar='THREADS', default='20')

    return expHTS_parser

def spadesParser(subparser):

    spades_parser = subparser.add_parser('spades', help='runs spades assembler')
    spades_parser.add_argument('-f', '--samplesfile', help='The filename of the sample file [default samples.txt]', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt')
    spades_parser.add_argument('-r', '--readFolder', help='Directory where the sequence data is stored [default 02-Cleaned]', action='store', dest='readFolder', default="02-Cleaned")
    spades_parser.add_argument('-n', '--spadesFolder', help='Directory where Spades output will be stored [defualt 03-Spades]', action='store', type=str, dest="spadesFolder", metavar='DIRECTORY', default='03-Spades')
    spades_parser.add_argument('-p', '--processors', help='Number of processors to be used by Spades [default 10]', action='store', type=str, dest="threads", metavar="PROCESSORS", default="10")
    spades_parser.add_argument('-e', '--error-correction', help='Perform error correction prior to assembly [default FALSE]', action='store_true', dest="errorCorrection", default=False)
    spades_parser.add_argument('-l', '--large_contigs_size', help='Size of contig considered "large" contigs [default 500]', action='store', type=str, dest="largeContigs", metavar="LARGE_CONTIG")
    spades_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default FALSE]', action='store_true', dest='force', default=False)
    spades_parser.add_argument('-m', '--memory', help='RAM limit for spades  [default: %(default)s]', action='store', type=str, dest="memory", metavar="RAM", default="250")
    
    return spades_parser

def fpParser(subparser):

    fp_parser = subparser.add_parser('forcepairs', help='creates new directory with forced pairs (split SE reads into R1 and R2)')
    fp_parser.add_argument('-f', '--samplesfile', help='The filename of the sample file [default: %(default)s]', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt')
    fp_parser.add_argument('-r', '--readFolder', help='Directory where the sequence data is stored [default: %(default)s]', action='store', dest='readFolder', default="02-Cleaned")
    fp_parser.add_argument('-n', '--forcepairFolder', help='Directory where Forced Pairs output will be stored [defualt: %(default)s]', action='store', type=str, dest="forcepairFolder", metavar='DIRECTORY', default='03-ForcedPairs')
    fp_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default: %(default)s', action='store_true', dest='force', default=False)
    
    return fp_parser

def kmerParser(subparser):

    kmer_parser = subparser.add_parser('kmernorm', help='Normilizes kmers')
    kmer_parser.add_argument('-f', '--samplesfile', help='The filename of the sample file [default: %(default)s]', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt')
    kmer_parser.add_argument('-r', '--readFolder', help='Directory where the sequence data is stored [default: %(default)s]', action='store', dest='readFolder', default="02-Cleaned")
    kmer_parser.add_argument('-k', '--kmernormFolder', help='Directory where Kmer output will be stored [default: %(default)s]', action='store', type=str, dest="kmerFolder", metavar='DIRECTORY', default='03-Kmer_Filter')
    kmer_parser.add_argument('-d', '--depth', help='Normilization depth [default: %(default)s]', action='store', type=str, dest="depth", metavar='DEPTH', default='20')
    kmer_parser.add_argument('-l', '--k_len', help='Kmer length [default: %(default)s]', action='store', type=str, dest="k_len", metavar='KMER LENGTH', default='15')
    kmer_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default: %(default)s]', action='store_true', dest='force', default=False)
    
    return kmer_parser

def parseArgs():
    revision_date = "Feb152016"
    parser = argparse.ArgumentParser(description="expHTS: Analysis of high throughput sequencing data in an experiment context ", epilog="For questions or comments, please contact Matt Settles <msettles@uidaho.edu>", add_help=True)
    parser.add_argument("--version", action="version", version="expHTS Version v" + version_num + "." + revision_date )
    subparsers = parser.add_subparsers(help='commands', dest='command')

    preprocessParser(subparsers)
    mappingParser(subparsers)
    htseqParser(subparsers)
    spadesParser(subparsers)
    fpParser(subparsers)
    kmerParser(subparsers)

    args = parser.parse_args()
    return args

def main():
    """
    main function
    """
    preprocess = preprocessCMD()
    mapping = mappingCMD()
    htseq = htseqCMD()
    spades = spadesCMD()
    forcepairs = forcepairCMD()
    kmernorm = kmerFilterCMD()
    commands = {'preprocess': preprocess, 'mapping': mapping, 'htseq': htseq, "spades": spades, "forcepairs": forcepairs, 'kmernorm': kmernorm}

    args = parseArgs()

    commands[args.command].execute(args)



signal.signal(signal.SIGINT, signal_handler)
