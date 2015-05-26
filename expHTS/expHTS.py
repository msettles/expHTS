# preprocess_app.py
import sys
import argparse
from bashSub import bashSub
from validate_app import validateApp
import signal
import os
from distutils import spawn


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
pgid = os.getpgid(os.getpid())

#the signal handler to kill the entire process group
#incase Cntr + C is sent (sub processes aren't killed
def signal_handler(signal, frame):
        print "Cntr + c was hit - ending process group number " + str(pgid)
        os.killpg(pgid, 9);



def preprocessParser(subparser):
	expHTS_parser = subparser.add_parser('preprocess', help='runs the baseline expHTS pipeline')
	expHTS_parser.add_argument('-f', '--file', help='The filename of the sample file [default samples.txt', action='store', type=str, dest='samplesFile', metavar='FILENAME', default='samples.txt');
	expHTS_parser.add_argument('-d', '--directory', help='Directory where the raw sequence data is stored [defualt 00-RawData]', action='store', type=str, dest='samplesDirectory', metavar='DIRECTORY', default='00-RawData');
	expHTS_parser.add_argument('-q', '--quality', help='Quality score to use during lucy trimming [default 20]', action='store', type=str, dest='qualityTrim', metavar='QUALITY', default='20');
	expHTS_parser.add_argument('-m', '--miniumumLength', help='Discard reads less than minimum length [default 50]', action='store', type=str, dest='minLength', metavar='MINIUMUMLENGTH', default='50');
	expHTS_parser.add_argument('-o', '--overlap', help='Overlap parameter for flash [default 700]', action='store', type=str, dest='overlapFlash', metavar='OVERLAP', default='700');
	expHTS_parser.add_argument('-O', '--skip-overlap', help='do not perform the overlapping using flash [default FALSE]', action='store', type=bool, dest='skipFlash', metavar='OVERLAP', default=False);
	expHTS_parser.add_argument('-s', '--skip-duplicates', help='do not preform the deduplication step [default FALSE]', action='store', type=bool, dest='skipDup', metavar='OVERLAP', default=False);
	expHTS_parser.add_argument('-c', '--contaminates-folder', help='folder name with contaminate sequences in fasta format [default NULL]', action='store', type=str, default="",  dest='contaminateFolder', metavar='CONTAMNANTS-FOLDER');
	expHTS_parser.add_argument('-F', '--final-folder', help='folder name in which the sequences will go [default 02-Cleaned', action='store', type=str, default="02-Cleaned", dest="finalDir", metavar='DIRECTORY')
	expHTS_parser.add_argument('-a', '--polyA', help='perform polyA trimming in sickle [default FALSE', action='store', type=bool, dest='polyTrim', default=False);
	#expHTS_parser.add_argument('-h', '--help', help='Show this help message and exit', action='store', type=bool, dest='helpOut', default=False);

	return expHTS_parser


def checkPreprocessApplications():
	applications = ["contaminant_screen.sh", "extract_unmapped_reads.py", "super_deduper", "sickle", "flash2"]
	
	for app in applications:
		if spawn.find_executable(app) is None:
			sys.stderr.write("It doesn't look like you have app - " + app + "\n" )
			exit(0)
		else:
			sys.stderr.write(app + " found\n")



def setupContaminateScreen(args):
	pass

def setupSuperDeduper(args):
	pass

def setupSickle2(args):
	pass

def setupFlash(args):
	pass

def returnReads(dictSampleSeqFiles):
	SE = ""
	PE1 = ""
	PE2 = ""

	#data struct 
	# { (sampleKey, seqKey) : [[SE], [SE], [PE1, PE2], [PE1, PE2]] }
	#diving into each of the sub lists in the dictionary value key
	for e in dictSampleSeqFiles:
		#if sublist only has one elment then it is SE read
		print e
		if len(e) == 1:
			if SE == "":
				SE = e[0]
			else:
				SE += "," + e[0]

		else:
			if PE1 == "":
				PE1 = e[0]
				PE2 = e[1]
			else:
				PE1 += "," + e[0]
				PE2 += "," + e[1]
				

	return [SE, PE1, PE2]


def check_dir(Dir):

	if not os.path.exists(Dir):
		os.mkdir(Dir)
		


class preprocessCMD:
    """
    preprocessing application for expHTS experiments
    """
    
    def __init__(self):
    	self.metaDataFolder = "MetaData"
        pass

    def execute(self, args):
        # ----------------------- options input files -----------------------
        print args
	checkPreprocessApplications()
	validate = validateApp()
	validate.setValidation(True)
	dictSampleSeqFiles = validate.validateSampleSheet(args.samplesDirectory, args.samplesFile)
	
	for key in dictSampleSeqFiles:
		check_dir(args.finalDir)
		check_dir(os.path.join(args.finalDir, key[1].rstrip()))
		check_dir(os.path.join(args.finalDir, self.metaDataFolder))
		
		meta = os.path.join(args.finalDir, self.metaDataFolder)

		SEandPE = returnReads(dictSampleSeqFiles[key])
			
		if SEandPE[0] != "":
			contArgsBaseline = "-e extract_unmapped_reads.py -o stdout -t 30 "
			if args.contaminateFolder != "":
				contArgsBaseline = "-c " + args.contaminateFolder +  contArgsBaseline
			cFilter = bashSub("contaminant_screen.sh", [SEandPE[0]], ['-U'], contArgsBaseline, os.path.join(meta, "SE_filter_info.log"))
			deduper = bashSub("super_deduper", cFilter.processSub(), ['-U'], "-p stdout", os.path.join(meta, "SE_deduper_info.log"))

			sickleArgs = " -o " + os.path.join(args.finalDir, key[1], "SE_not_merged.fastq")  +  " -t sanger -l " + args.minLength
			if args.polyTrim:
				sickleArgs += " -a "
			
			sickle =  bashSub("sickle se", deduper.processSub(), ['-f'], sickleArgs, "SE_sickle_info.log")
			print "___ SE COMMANDS ____"
			print sickle.getCommand()
			sickle.runCmd("")

		if SEandPE[1] != "":
			contArgsBaseline = "-e extract_unmapped_reads.py -o stdout -t 30 "
			if args.contaminateFolder != "":
				contArgsBaseline = "-c " + args.contaminateFolder +  contArgsBaseline
			cFilter = bashSub("contaminant_screen.sh", [SEandPE[1], SEandPE[2]], ['-1', '-2'], contArgsBaseline, os.path.join(meta, "PE_filter_info.log"))
			deduper = bashSub("super_deduper", cFilter.processSub(), ['-i'], "-p stdout", os.path.join(meta, "PE_deduper_info.log"))

			sickleArgs = " -M stdout -t sanger -l " + args.minLength
			if args.polyTrim:
				sickleArgs += " -a "
	
			sickle =  bashSub("sickle pe", deduper.processSub(), ['-c'], sickleArgs , os.path.join(meta, "PE_sickle_info.log"))
			flash = bashSub("flash2", sickle.processSub(), ['--interleaved-input'], " --allow-outies -d " + os.path.join(args.finalDir, key[1]) + " ", os.path.join(meta, "flash_info.log"))
			print "___ PE COMMANDS ___"
			print flash.getCommand()
			flash.runCmd("")

        return 0


	

def parseArgs():

	parser = argparse.ArgumentParser(description = "expHTS is a python application that is awesome", epilog="For questions or comments, please contact Matt Settles <msettles@uidaho.edu>", add_help=True)
	parser.add_argument("--version", action="version", version="%(progs)s Version " + version_num)
	subparsers = parser.add_subparsers(help='commands', dest='command')
	preprocessParser(subparsers)
	args = parser.parse_args()
	return args


def main():
	parser = argparse.ArgumentParser(description='Runs expHTS pipline', epilog='Needs update', add_help=True)
	subparsers = parser.add_subparsers(help='commands', dest='command')
	preprocess = preprocessCMD()
	commands = {'preprocess': preprocess}
	args = parseArgs()
	commands[args.command].execute(args)
	
	#return args

signal.signal(signal.SIGINT, signal_handler)
main()
