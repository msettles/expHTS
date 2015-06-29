#!/usr/bin/python


from subprocess import PIPE, STDOUT, Popen, check_call
import os
import time
import sys
import signal
import optparse
import time

flash_path = "/mnt/home/grcuser/module_grc/src/ETS_tmp/flash2" 
sickle_path = "/mnt/home/grcuser/module_grc/src/ETS_tmp/sickle"
super_deduper_path = "/mnt/home/grcuser/module_grc/src/ETS_tmp/super_deduper"
screen_path = "/mnt/home/grcuser/module_grc/src/ETS_tmp"


#Basic Application class
class bashSub:
	app = ""
	files = ""
	file_args = ""	
	args = ""
	f_info = ""	
	cmd = ""

	def __init__(self, app, files, file_args, args, info_file):
		self.app = app
		self.files = files
		self.file_args = file_args
		self.args = args
		self.f_info = info_file
		self.time = 0

		if len(files) != len(file_args):
			print >> sys.stderr, "Error: File args need to be the same length"

		self.cmd = app
		
		for i in range(0, len(files)):
			self.cmd += " " + file_args[i] + " " + files[i] + " "
		
		self.cmd += " " + args

		if self.f_info != "":
			self.cmd += " 2>" + self.f_info
		
	
	def getCommand(self):
		return self.cmd

	def runCmd(self, additional_args):
		try:
			start = time.time()
			p2 = Popen(self.cmd + additional_args, stdout=sys.stdout, stdin=PIPE, stderr=sys.stderr, shell=True, executable = "/bin/bash");
			(output, error) = p2.communicate()
			sys.stderr.flush();
			end = time.time()
			self.time = end - start
			print "Seconds: " + str(self.time)
		except:
			print "error"
			print self.cmd + additional_args
			print "Command failed"

	def returnTime(self):
		return self.time


	def processSub(self):
		return ["<(" + self.cmd + ")"]




def main():
	
	parser = optparse.OptionParser()
	
	parser.add_option("-1", "--read1", dest="r1_file", default="")
	parser.add_option("-2", "--read2", dest="r2_file", default="")
	parser.add_option("-U", "--single-end", dest="single_end", default="")
	parser.add_option("-I", "--interleaved", dest="interleaved", default="")
	parser.add_option("-a", "--args", dest="args", default="")
	parser.add_option("-o", "--output", dest="output", default="")

	(options, args) = parser.parse_args()
	output = os.path.dirname(options.output)


	single = False
	if options.r1_file != "" and options.r2_file != "":
		reads = [options.r1_file, options.r2_file]
		file_args = ['-1', '-2']
		sickle_app = " pe "
		sickle_out = " -M"
		sickle_file = ['-c']
		deduper_file = ['-i']
		print "Read 1 and Read 2"
	elif options.single_end != "":
		single = True
		reads = [options.single_end]
		file_args = ['-U']
		sickle_app = " se "
		sickle_file = ['-f']
		sickle_out = " -o"
		deduper_file = ['-U']
		print "Single end"
	elif options.interleaved != "":
		reads = [options.interleaved]
		file_args = ['-I']
		sickle_file = ['-c']
		sickle_args = " -o "
		print "Interleaved"
	else:
		print "No interleaved, single end, nor paired end files specified. Exiting."
		exit(0)

	print "Starting up the processes"
	#print "Starting Contaminate Screening\nstderr going to " + "screen_info.txt"
	c_filter = Application(screen_path, reads, files_args, options.args + " -e /mnt/home/stre3949/ETS_Pipeline/extract_unmapped_reads.py -o stdout", output + "screen_info.txt")
	print "Starting Super Deduper\nstderr going to " + "deduper_info.txt"
	deduper = Application(super_deduper_path, reads, file_args, "-p stdout", output + "deduper_info.txt");
	print "Starting Sickle\nstderr going to " + output + "sickle_info.txt"
	sickle = Application(sickle_path + sickle_app, deduper.processSub(), sickle_file, sickle_out + " stdout -t sanger -l 150 ", output + "sickle_info.txt")
	#sickle = Application(sickle_path + sickle_app, c_filter.processSub(), sickle_file, sickle_out + " stdout -t sanger -a -l 35 ", output + "sickle_info.txt")

	if not single:
		print "Starting Flash\nstderr going to " + output + "flash.txt"
		print "Final Output out.extendedFrags.fastq, out.notComined_1.fastq, out.notCombined_2.fastq"
		flash = Application(flash_path, sickle.processSub(), ["--interleaved-input"], " -d " + output,  output + "flash.txt")
		flash.runCmd("")
		print flash.getCommand()
	else:
		print "Final Output " + output
		sickle.runCmd(">" + output) 

	print "Everything is finished up now"

