import re
import os

def flashParser(f, h, d):
	if not os.path.isfile(f):
		return

	f = open(f, "r")
	strSearch = ["[FLASH]     Total pairs:", "[FLASH]     Discarded pairs:", "[FLASH]     Combined pairs:", "[FLASH]         Innie pairs:", "[FLASH]         Outie pairs:", "[FLASH]     Uncombined pairs:", "[FLASH]     Percent combined:"]

	header = ["Flash_Pairs", "Flash_Discarded", "Combined_Pairs", "Innie_Combined", "Innie_Combined_Percent",  "Outie_Combined", "Outie_Combined_Percent", "Uncombined", "Percent_Combined"]

	data = []
	lines = f.readlines()
	if lines == []:
		return;

	index = 0
	for e in lines:
		if strSearch[index] in e:
			if strSearch[index] == "[FLASH]     Percent combined:":
				data.append(float(re.findall("\d+.\d+", e)[-1]));
			else:
				data.append(int(re.search(r'\d+', e).group()))

			if "[FLASH]         Innie pairs:" == strSearch[index] or "[FLASH]         Outie pairs:" == strSearch[index]:
				data.append(float(re.findall("\d+.\d+", e)[-1]));

			index += 1
			if index >= len(strSearch):
				break;
	h += header
	d += data
def sickleParser(f, h, d):
	if not os.path.isfile(f):
		return

	f = open(f, "r")
	
	strSearch = ["FastQ paired records kept:", "FastQ single records kept:", "FastQ paired records discarded:", "FastQ single records discarded:", "PE1 Base pairs left removed:", 
		"PE1 Base pairs right removed:", "PE2 Base pairs left removed", "PE2 Base pairs right removed:"]

	header  = ["Sickle_Pairs_Kept", "Sickle_Single_Kept", "Sickle_Pairs_Discarded", "Sickle_Single_Discarded", "Sickle_R1_BP_Left_Removed", "Sickle_R1_BP_Right_Removed", "Sickle_R2_BP_Left_Removed", "Sickle_R2_BP_Right_Removed"]
	headerPolyA = ["Sickle_R1_PolyAT_Trim", "Sickle_R1_PolyAT_Trim"]

	data = []
	dataPolyA = []

	lines = f.readlines()
	index = 0
	if lines == []:
		return;
	
	for e in lines:
		if strSearch[index] in e:
			e = e.split(":")[1]
			data.append(int(re.search(r'\d+', e).group()))
			index += 1
			if index >= len(strSearch):
				break;
		if "Poly AT tails" in e:
			e = e.split(":")[1]
			dataPolyA.append(int(re.search(r'\d+', e).group()))
			
		
	h += header
	d += data
	
	if len(dataPolyA) != 0:
		h += headerPolyA
		d += dataPolyA		


def scythParser(f, h, d):
	if not os.path.isfile(f):
		return
	f = open(f, "r")

	strSearch = ["contaminated: "]
	header = ["Scythe_contaminated", "Scythe_no_contamination", "Scythe_total"]
	data = []

	lines = f.readlines()
	if lines == []:
		return;

	index = 0;
	
	for e in lines:
		if strSearch[index] in e:
			info = e.split(',')
			data.append(int(re.search(r'\d+', info[0]).group()))
			data.append(int(re.search(r'\d+', info[1]).group()))
			data.append(int(re.search(r'\d+', info[2]).group()))
			break;
	
	h += header
	d += data

def filterParser(f, h, d):
	if not os.path.isfile(f):
		return
	f = open(f, "r")
		

	strSearch = ["PE_written:", "SE_written:", "Discarded:",]
	
	header = ["Screeing_PE_writtened", "Screening_SE_written", "Screening_discarded"]
	data = []

	lines = f.readlines()
	index = 0
	if lines == []:
		return;
	
	lines = lines[-1]
	lines = lines.split("|")
	
	for e in lines:
		if strSearch[index] in e:
			data.append(int(re.search(r'\d+', e).group()))
			index += 1
			if index >= len(strSearch):
				break;
			
	

	h += header
	d += data



def deduperParser(f, h, d):
	if not os.path.isfile(f):
		return
	f = open(f, "r")
		

	strSearch = ["reads:", "duplicates:", "percent:", "discarded:"]
	
	header = ["Deduper_reads", "Deduper_duplicates", "Deduper_Percent_Duplicate", "Deduper_discarded"]
	data = []

	lines = f.readlines()

	if lines == []:
		return;

	lines = lines[-1]
	lines = lines.split("|")
	index = 0

	for e in lines:
		if strSearch[index] in e:
			if strSearch[index] == "percent:":
				data.append(float(re.findall("\d+.\d+", e)[-1]))
			else:
				data.append(int(re.search(r'\d+', e).group()))
			index += 1
			if index >= len(strSearch):
				break;
			
	
	h += header
	d += data

def finalParser(f, h, d):
	if not os.path.isfile(f):
		return
	f = open(f, "r")
		

	lines = f.readlines()

	if lines == []:
		return;
	
	h += lines[0].split("\t")	
	d += lines[1].split("\t")	


def printToFile(out, header, data):
	f = open(out, "w")
	
	f.write('\t'.join(map(str,header)))
	f.write('\t'.join(map(str,data)))

	f.close()

def parseOut(base, sample):
	import os
	data = []
	header = []

	
	header.append("Sample")
	data.append(sample)


	filterParser(os.path.join(base,"PE_filter_info.log"), header, data)
	deduperParser(os.path.join(base, "PE_deduper_info.log"), header, data)
	sickleParser(os.path.join(base, "PE_sickle_info.log"), header, data)
	flashParser(os.path.join(base, "flash_info.log"), header, data)
	finalParser(os.path.join(base, "finalCleanup.log"), header, data)

	printToFile(os.path.join(base, sample + "_SummaryStats.log"), header, data)	

	return os.path.join(base, sample + "_SummaryStats.log")

def bringTogether(listFiles, out):
	first = 0
	outFile = open(out, "w")

	for e in listFiles:
	
		f = open(e, "r")
		lines = f.readlines()	
		if first == 0:
			outFile.write(lines[0])
			outFile.write(lines[1])
		else:
			outFile.write(lines[1])
		first += 1	
		f.close()

	outFile.close()
