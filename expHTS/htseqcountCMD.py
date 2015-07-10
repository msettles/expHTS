from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOut, bringTogether
from bashSub import bashSub

def checkPreprocessApplications():
	applications = ["./contaminant_screen.sh", "./extract_unmapped_reads.py", "super_deduper", "sickle", "flash2"]

	for app in applications:
		if spawn.find_executable(app) is None:
			sys.stderr.write("It doesn't look like you have app - " + app + "\n" )
			exit(0)
		else:
			sys.stderr.write(app + " found\n")



def returnReads(dictSampleSeqFiles):
	SE = ""
	PE1 = ""
	PE2 = ""

#data struct 
# { (sampleKey, seqKey) : [[SE], [SE], [PE1, PE2], [PE1, PE2]] }
#diving into each of the sub lists in the dictionary value key
	for e in dictSampleSeqFiles:
		#if sublist only has one elment then it is SE read
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



class htseqCMD:
	def __init__(self):
		self.metaDataFolder = "MetaData"

		
	def index(self, ref):
		if not os.path.exists(ref):
			print "Would you mind adding a gtf file? (-R) Thank you."
			exit(1);


	def execute(self, args):
		logFiles = []
		checkPreprocessApplications()
		validate = validateApp()
		validate.setValidation(True)
		self.index(args.refGTF)
		dictSampleSeqFiles = validate.validateSampleSheetHTSeq(args.readFolder, args.finalDir, args.samplesFile, args.force)
		for keys in dictSampleSeqFiles.keys():
			check_dir(args.finalDir)
			check_dir(keys[1])

			bamFile = os.path.join(keys[0], keys[0].split("/")[-1]) + ".bam"
			outFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".out"
			countFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".counts"

			cmdString = "htseq-count -f bam -r " + args.order + " -s " + args.stranded + " -a 10 -t exon -i gene_id -m " + args.mode + " " + bamFile + " " + args.refGTF + " 2>" + outFile + " >" + countFile

			htseqCmd = bashSub(cmdString, [''], [''], '', '')
			print htseqCmd.getCommand()
			htseqCmd.runCmd("")
			
		time = 0
		print dictSampleSeqFiles

        def clean(self):
                import glob
                for f in glob.glob(".screening_cont*"):
                        os.remove(f)


