# preprocess_app.py
import sys
import argparse
from bashSub import bashSub
from parse_files import parseOut, bringTogether
import signal
import os
from distutils import spawn
from preprocessCMD import preprocessCMD
from mappingCMD import mappingCMD

version_num = "0.0"
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

#get process group number to know what to kill
#get process group number to know what to kill
pgid = os.getpgid(os.getpid())

#the signal handler to kill the entire process group
#incase Cntr + C is sent (sub processes aren't killed
def signal_handler(signal, frame):
        print "Cntr + c was hit - ending process group number " + str(pgid)
        import glob
        for f in glob.glob(".screening_cont*"):
                os.remove(f)
        os.killpg(pgid, 9);


def mappingParser(subparser):

        mapping_parser = subparser.add_parser('mapping', help='runs the baseline expHTS pipeline')
        mapping_parser.add_argument('-f', '--file', help='The filename of the sample file [default samples.txt', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt');
        mapping_parser.add_argument('-r', '--readFolder', help='Directory where the sequence data is stored [defualt 02-Cleaned]', action='store', type=str, dest='readFolder', metavar='FOLDER', default='02-Cleaned');
        mapping_parser.add_argument('-R', "--reference", help='Reference fasta to map against', action='store', type=str, dest='refFasta', metavar='REFERENCE', default='');
        mapping_parser.add_argument('-M', "--mappingAlgorithm", help='Mapping algorithm bwa or bowtie2 [defualt bwa]', action='store', type=str, dest='mapping', metavar='ALGORITHM', default='bwa');
        mapping_parser.add_argument('-t', '--threads', help='Number of threads to be used [Default 20]', action='store', type=str, dest='threads', metavar='THREADS', default='20');
        mapping_parser.add_argument('-F', '--final-folder', help='folder name in which the sequences will go [default 03-Cleaned]', action='store', type=str, default="03-BWA", dest="finalDir", metavar='DIRECTORY')
        mapping_parser.add_argument('-a', '--polyA', help='perform polyA trimming in sickle [default FALSE]', action='store_true', dest='polyTrim', default=False)
	mapping_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default FALSE]', action='store_true', dest='force', default=False)
	mapping_parser.add_argument('-i', '--force-index', help='overwrites old index files [default FALSE]', action='store_true', dest='forceIndex', default=False)

        return mapping_parser


def preprocessParser(subparser):

        expHTS_parser = subparser.add_parser('preprocess', help='runs the baseline expHTS pipeline')
        expHTS_parser.add_argument('-f', '--file', help='The filename of the sample file [default samples.txt', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt');
        expHTS_parser.add_argument('-S', '--forceSplit', help='Forces splits of SE reads [default FALSE]', action='store_true', dest='forceSplit', default=False);
        expHTS_parser.add_argument('-t', '--threads', help='Threads for bowtie2 [default 20]', action='store', type=str, dest='threads', metavar='THREADS', default='20');
        expHTS_parser.add_argument('-A', '--adapterfasta', help='folder name with adapter sequences in fasta format [default truseq adapter sequence]', action='store', type=str, default=r'<(printf ">TruSeq_forward_contam\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[NNNNNN]ATCTCGTATGCCGTCTTCTGCTTGAAAAA\n>TruSeq_reverse_contam\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA")',  dest='adapter', metavar='CONTAMNANTS-FOLDER');
        expHTS_parser.add_argument('-d', '--directory', help='Directory where the raw sequence data is stored [defualt 00-RawData]', action='store', type=str, dest='samplesDirectory', metavar='DIRECTORY', default='00-RawData');
        expHTS_parser.add_argument('-q', '--quality', help='Quality score to use during lucy trimming [default 20]', action='store', type=str, dest='qualityTrim', metavar='QUALITY', default='20');
        expHTS_parser.add_argument('-m', '--miniumumLength', help='Discard reads less than minimum length [default 50]', action='store', type=str, dest='minLength', metavar='MINIUMUMLENGTH', default='50');
        expHTS_parser.add_argument('-o', '--overlap', help='Overlap parameter for flash [default 700]', action='store', type=str, dest='overlapFlash', metavar='OVERLAP', default='700');
        expHTS_parser.add_argument('-O', '--skip-overlap', help='do not perform the overlapping using flash [default FALSE]', action='store_true', dest='skipFlash',  default=False);
        expHTS_parser.add_argument('-s', '--skip-duplicates', help='do not preform the deduplication step [default FALSE]', action='store_true',  dest='skipDup',  default=False);
        expHTS_parser.add_argument('-c', '--contaminates-folder', help='folder name with contaminate sequences in fasta format [default NULL]', action='store', type=str, default='',  dest='contaminateFolder', metavar='CONTAMNANTS-FOLDER');
        expHTS_parser.add_argument('-F', '--final-folder', help='folder name in which the sequences will go [default 02-Cleaned', action='store', type=str, default="02-Cleaned", dest="finalDir", metavar='DIRECTORY')
        expHTS_parser.add_argument('-a', '--polyA', help='perform polyA trimming in sickle [default FALSE]', action='store_true', dest='polyTrim', default=False);
        expHTS_parser.add_argument('-w', '--overwrite', help='overwrite a sequence id folder [default FALSE]', action='store_true', dest='force', default=False);

        return expHTS_parser


def parseArgs():
        parser = argparse.ArgumentParser(description = "expHTS is a python application that is awesome", epilog="For questions or comments, please contact Matt Settles <msettles@uidaho.edu>", add_help=True)
        parser.add_argument("--version", action="version", version="%(progs)s Version " + version_num)
        subparsers = parser.add_subparsers(help='commands', dest='command')
        preprocessParser(subparsers)
        mappingParser(subparsers)
        args = parser.parse_args()
        return args

def main():
        parser = argparse.ArgumentParser(description='Runs expHTS pipline', epilog='Needs update', add_help=True)
        subparsers = parser.add_subparsers(help='commands', dest='command')
        preprocess= preprocessCMD()
        mapping = mappingCMD()
	print "Here"
        commands = {'preprocess': preprocess, 'mapping': mapping}
        args = parseArgs()
        print args.command
        commands[args.command].execute(args)

signal.signal(signal.SIGINT, signal_handler)

