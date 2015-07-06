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



class mappingCMD:
	def __init__(self):
		self.metaDataFolder = "MetaData"

	def createIndex(self, ref, algorithm):
		if os.path.exists(ref):
			if "bowtie" in algorithm:
				if not os.path.exists(os.path.join(ref, ".bt2")):
					createTheIndex = bashSub("bowtie2-build ", [ref], [''], '', '/dev/null')
					print createTheIndex.getCommand()
					createTheIndex.runCmd("");
			else:
				if not os.path.exists(os.path.join(ref, ".sa")):
					createTheIndex = bashSub("bwa index ", [ref], [''], '', '/dev/null')
					print createTheIndex.getCommand()
					createTheIndex.runCmd("");

		else:
			print "Doesn't seem ref - " + ref + "actually exists"
		
		
	def index(self, ref, algorithm):
		if ref == "":
			print "Would you mind adding a reference file? (-R) Thank you."
		else:
			self.createIndex(ref, algorithm)
			

	def execute(self, args):
		logFiles = []
		checkPreprocessApplications()
		validate = validateApp()
		validate.setValidation(True)
		dictSampleSeqFiles = validate.validateSampleSheet(args.readFolder, args.finalDir, args.samplesFile, args.force, True)
		time = 0

		self.index(args.refFasta, args.mapping)

		for key in dictSampleSeqFiles:
			check_dir(args.finalDir)
			check_dir(key[1])
			meta =  key[1]

			endString = '| tee >(grep "^@" >headers.log) | tee >(samtools flagstat - >flagstats.log) | samtools view -bS - >' + os.path.join(key[1], "test.bam")
                        SEandPE = returnReads(dictSampleSeqFiles[key])

			if SEandPE[0] != "":
				terminalString = []
			if SEandPE[1] != "":
				terminalString = []

				terminalString.append(bashSub("bwa mem", [args.threads], ['-t'], args.refFasta + " " + SEandPE[1] + " " + SEandPE[2] + " -M " + endString, "/dev/null"))
				
				print "___ PE COMMANDS ___"
				print terminalString[-1].getCommand()
				terminalString[-1].runCmd("")
				sys.stderr.flush()
				time += terminalString[-1].returnTime()
				logFiles.append(parseOut(key[1], key[1].split("/")[-1]))

	                #bringTogether(logFiles, os.path.join(key[1].split("/")[0], "stats.log"))
       	        	print "Total amount of seconds to run all samples"
                	print "Seconds: " + str(time)

                	self.clean()


        def clean(self):
                import glob
                for f in glob.glob(".screening_cont*"):
                        os.remove(f)


