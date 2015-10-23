from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOutHTseq, bringTogether
from bashSub import bashSub


def checkPreprocessApplications():
    applications = ["spades.py"]
    source = ["http://bioinf.spbau.ru/spades"]
    i = 0;
    for app in applications:
        if spawn.find_executable(app) is None:
            sys.stderr.write("It doesn't look like you have app - " + app + "\n")
            sys.stderr.write("Download it here - " + source[i] + "\n");
            exit(0)
        else:
            sys.stderr.write(app + " found\n")
            i += 0

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


class spadesCMD:

    def __init__(self):
        self.metaDataFolder = "MetaData"


    def execute(self, args):
        time =  0
        checkPreprocessApplications();
        logFiles = []

        # checkPreprocessApplications()
        validate = validateApp()
        validate.setValidation(True)
        dictSampleSeqFiles = validate.validateSampleSheet(args.readFolder, args.spadesFolder, args.samplesFile, args.force, True)

        for keys in dictSampleSeqFiles.keys():
            check_dir(args.spadesFolder)
            check_dir(keys[1])
            terminal = []
            #countFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".counts"

            if (len(dictSampleSeqFiles[keys][0]) == 3):
                terminal.append(bashSub("spades.py", dictSampleSeqFiles[keys][0], ['-1', '-2', '-s'], " --careful -t " + args.threads + " -o " + args.spadesFolder + " -m " + args.memory, ''))
            elif (len(dictSampleSeqFiles[keys][0]) == 2):
                terminal.append(bashSub("spades.py", dictSampleSeqFiles[keys][0], ['-1', '-2'], "--careful -t " + args.threads + " -o " + args.spadesFolder + " -m " + args.memory, ''))

            print terminal[-1].getCommand()
            terminal[-1].runCmd("")
            sys.stderr.flush()
            #time += runSortByName.returnTime() + runView.returnTime() + htseqCmd.returnTime()

            #logFiles.append(parseOutHTseq(keys[1], keys[1].split("/")[-1]))

        #bringTogether(logFiles, os.path.join(args.finalDir, "Counts_Summary.log"))

        print "Total amount of seconds to run all samples"
        print "Seconds: " + str(time)

