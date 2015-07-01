import sys
import os

class validateApp:
	
	sampleFiles = {}

	def __init__(self):
		self.exitValidate = False
		self.debug = True
		self.verbose = True


	def validateSampleSheet(self, dirSample, finalDir, sampleSheet, force):
		linenum = 0
		if not os.path.exists(dirSample):
			self.exitTime("Directory " + dirSample + " is not found")
		if os.path.exists(sampleSheet):
			f = open(sampleSheet, "r");
		else:
			self.exitTime("No sample sheet named " + sampleSheet + " was found")

		
		for e in f.readlines():
			if linenum != 0:
				self.sampleSequenceID(dirSample, finalDir, e.split("\t"), force)
			linenum += 1

		for key in self.sampleFiles:
			self.sampleFiles[key].sort(key=lambda x: len(x))

		return self.sampleFiles

	def sampleSequenceID(self, dirSample, finalDir, seqID, force):
		if seqID[0][0] == "#":
			pass
		elif len(seqID) >= 2:

			seqID[0] = os.path.join(dirSample, seqID[0])
			seqID[1] = os.path.join(finalDir, seqID[1].rstrip())

			self.finalDirTest(seqID[1], force)
			self.directoryFiles(seqID)
		else:
			self.exitTime("There wasn't two columns in the file")

	def finalDirTest(self, sampleID, force):
		if os.path.exists(sampleID) and not force:
			self.exitTime(sampleID + " was all ready created. Use the -w or --overwrite option to overwrite")
		elif os.path.exists(sampleID):
			print "Warning"
			print "Overwrite was turned on - overwriting " + sampleID + "\n\n"
			
			

	#sets up the directory dictionary
	#set up key with tuple (sample and seq)
	def directoryFiles(self, sampleSeq):
		sampleSeq = tuple(sampleSeq)
		#print sampleSeq
		directoryTest = sampleSeq[0]
		#fastqCount insures at least one fastq files is under the directory
		fastqCount = 0

		if self.testDirectory(directoryTest):
			for subdir, dir, files in os.walk(directoryTest):
				for file in files:
					file = os.path.abspath(os.path.join(directoryTest, file))
					if "_R1" in file and ".fastq" in file:
						fastqCount += 1
						if not sampleSeq in self.sampleFiles:
							self.sampleFiles[sampleSeq] = []
						
						self.sampleFiles[sampleSeq].append(self.isPairedEnd(file))


			if (fastqCount == 0):
				self.exitTime("No fastq files were found under this directory - " + directoryTest)

		else:
			self.exitTime("Direcotry " + directoryTest + " does not exists")
			
	def exitTime(self, exitString):
		print exitString
		if self.exitValidate:
			print "Failure in Validation - exiting process"
			exit(-1)
		else:
			print "Failure in Validation - going to continue to continue to validate sample sheet"

	#returns 2 elements in list if paired end reads or 1 element in list if it is a SE
	def isPairedEnd(self, fileRead1):
		file2 = fileRead1.replace("_R1", "_R2")
		if os.path.exists(file2):
			return [fileRead1, file2]
		else:
			return [fileRead1]
		

	#test if directory exists . . . meh. probably can get ride of it
	def testDirectory(self, directoryTest):
		return os.path.isdir(directoryTest)

	#outputs the single end paired end reads	
	def infoOutput(self):
		#print self.sampleFiles
		for key in self.sampleFiles:
			#print "SEQUENCE " + key[0] + " SAMPLE " + key[1]
			for files in self.sampleFiles[key]:
				if len(files) == 1:
					pass
					#print "SE Files: " + str(files).strip("[]")
				else:
					pass
					#print "PE Files: " + str(files).strip("[]")

	def setValidation(self, exitOut):
		self.exitValidate = exitOut

		
	#tuple key (seq ID - sample ID)
	#value - 1 value is SE reads and 2 values is PE
	def dictionaryFilesReturn(self):
		return self.sampleFiles

def main():
	test = validateApp(False)
	test.validateSampleSheet("samples.txt")
	test.infoOutput()
