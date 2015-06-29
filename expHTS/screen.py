from bashSub import bashSub
import os
import glob
import argparse


class screening:
	def __init__(self):
		pass


	def screen(self, reads, threads, screenFolder):
		#wget phix
		#concat everything in folder
		#bowtie2-build

		fastaFile = ".screening_cont.fasta"
		indexFile = ".screening_cont_index"

		wget = bashSub("wget", ["-"], ["-O"], '"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=9626372"', "/dev/null")
		cat = None

		if screenFolder != "":
			if os.path.exists(screenFolder):
				cat = bashSub("cat", [""], [""], os.path.join(screenFolder, "*.fa*"), "")
			else:
				print screenFolder + " does not exists! Please enter a proper folder"
				exit(0)



		if cat is None:
			combine = bashSub("cat", wget.processSub(), [""], "1>" + fastaFile, "/dev/null")
		else:
			combine = bashSub("cat", [wget.processSub(), cat.processSub()], ["", ""], "1>" + fastaFile, "/dev/null")
		


		if len(reads) == 3:
			files = ['-1', '-2', '-U']
		elif len(reads) == 2:
			files = ['-1', '-2']
		else:
			files = ['-U']

		build = bashSub("bowtie2-build", [fastaFile, indexFile], ["", ""], " 1>/dev/null ", "/dev/null")
		bowtie = bashSub("bowtie2", reads, files, "-x " + indexFile + " -X 1500 -I 0 --very-sensitive-local -p " + threads + " -q", "/dev/null")
		source = bashSub("source", ["/usr/modules/init/bash"], [""], "", "/dev/null")
		load = bashSub("module", [""], [""], "load bowtie2 grc/2.0", "/dev/null")

		source.runCmd("")
		load.runCmd("")

		if not os.path.exists(fastaFile):
			combine.runCmd("")
			build.runCmd("")

		cmd = bowtie.runCmd("")
	


	

	def cleanUp(self):
		for f in glob.glob(".screening_cont*"):
			os.remove(f)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-1", "--read1", type=str, action= "store", dest="read1",help="Read 1", default = "")
	parser.add_argument("-2", "--read2", type=str, action="store", dest="read2", help="Read 2", default = "")
	parser.add_argument("-U", "--single", type=str, action="store", dest="single", help="Single End", default = "")
	parser.add_argument("-t", "--threads", type=str, action="store", dest="threads", default="20", help="threads")
	parser.add_argument("-c", "--contimantFolder", type=str, action="store", dest="cont", help="location of fasta files you would like to remove", metavar="DIR", default="")
	args = parser.parse_args()	
	screen = screening()

	if (args.read1 == "" and args.read2 != "") or (args.read1 != "" and args.read2 == ""):
		print "Read 1 and Read 2 both must be specified if one is"
		exit(-1)
	elif args.read1 != "" and args.single != "":
		screen.screen([args.read1, args.read2, args.single], args.threads, args.cont)
	elif args.single != "":
		screen.screen([args.single], args.threads, args.cont)
	elif args.read1 != "":
		screen.screen([args.read1, args.read2], args.threads, args.cont)
	else:
		print "Please specifiy reads"
		exit(1)





main()
