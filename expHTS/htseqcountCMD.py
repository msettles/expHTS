from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOutHTseq, bringTogether
from bashSub import bashSub


def checkPreprocessApplications():
    applications = ["./contaminant_screen.sh", "./extract_unmapped_reads.py", "super_deduper", "sickle", "flash2"]

    for app in applications:
        if spawn.find_executable(app) is None:
            sys.stderr.write("It doesn't look like you have app - " + app + "\n")
            exit(0)
        else:
            sys.stderr.write(app + " found\n")


def returnReads(dictSampleSeqFiles):
    SE = ""
    PE1 = ""
    PE2 = ""

    # data struct
    # { (sampleKey, seqKey) : [[SE], [SE], [PE1, PE2], [PE1, PE2]] }
    # diving into each of the sub lists in the dictionary value key
    for e in dictSampleSeqFiles:
        # if sublist only has one elment then it is SE read
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
            exit(1)

    def execute(self, args):
        time =  0
        logFiles = []

        # checkPreprocessApplications()
        validate = validateApp()
        validate.setValidation(True)
        dictSampleSeqFiles = validate.validateSampleSheetHTSeq(args.readFolder, args.finalDir, args.samplesFile, args.force)

        self.index(args.refGTF)

        for keys in dictSampleSeqFiles.keys():
            check_dir(args.finalDir)
            check_dir(keys[1])

            bamFile = os.path.join(keys[0], keys[0].split("/")[-1]) + ".bam"
            outFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".out"
            countFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".counts"

            # runSortByName = bashSub("samtools view -bF 0x100", [bamFile], [''], "| samtools sort -n - " + os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid"), '/dev/null')
            runSortByName = bashSub("samtools sort -n", [bamFile], [''],  os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid"), '/dev/null')
            print runSortByName.getCommand()
            runSortByName.runCmd("")

            runSortByName = bashSub("samtools sort -n ", [os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid")], [''], os.path.join(keys[1], keys[1].split('/')[-1]) + ".byreadid" , '/dev/null')
            print runSortByName.getCommand()
            runSortByName.runCmd("")

            runView = bashSub("samtools view -F 0x100 ", [os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid.bam")], [''], "> " + os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid.sam"), '/dev/null')
            print runView.getCommand()
            runView.runCmd("")

            cmdString = "htseq-count -f sam -s " + args.stranded + " -a " + args.minQual + " -t " + args.type + " -i " + args.idattr + " -m " + args.mode + " " + os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid.sam ") + args.refGTF + " 2>" + outFile + " >" + countFile

            htseqCmd = bashSub(cmdString, [''], [''], '', '')
            print htseqCmd.getCommand()
            htseqCmd.runCmd("")

            sys.stderr.flush()
            time += runSortByName.returnTime() + runView.returnTime() + htseqCmd.returnTime()

            logFiles.append(parseOutHTseq(keys[1], keys[1].split("/")[-1]))

        bringTogether(logFiles, os.path.join(args.finalDir, "Counts_Summary.log"))

        print "Total amount of seconds to run all samples"
        print "Seconds: " + str(time)

