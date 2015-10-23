from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOutHTseq, bringTogether
from bashSub import bashSub


def checkPreprocessApplications():
    applications = ["kmer_filter"]
    source = ["http://catchenlab.life.illinois.edu/stacks/"]
    i = 0

    for app in applications:
        if spawn.find_executable(app) is None:
            sys.stderr.write("It doesn't look like you have app - " + app + "\n")
            sys.stderr.write("You can download it here - " +  source[i] + "\n")
            exit(0)
        else:
            sys.stderr.write(app + " found\n")
            i += 1

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


class kmerFilterCMD:

    def __init__(self):
        self.metaDataFolder = "MetaData"

    def index(self, ref):
        if not os.path.exists(ref):
            print "Would you mind adding a gtf file? (-R) Thank you."
            exit(1)

    def execute(self, args):
        time =  0
        checkPreprocessApplications();
        logFiles = []

        # checkPreprocessApplications()
        validate = validateApp()
        validate.setValidation(True)
        dictSampleSeqFiles = validate.validateSampleSheet(args.readFolder, args.kmerFolder, args.samplesFile, args.force, True)

        for keys in dictSampleSeqFiles.keys():
            check_dir(args.kmerFolder)
            check_dir(keys[1])
            terminal = []
            #countFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".counts"

            if (len(dictSampleSeqFiles[keys][0]) == 3):
                terminal.append(bashSub("kmer_filter", dictSampleSeqFiles[keys][0], ['-1', '-2', '-f'], "--normalize " + args.depth + " --k_len " + args.k_len + ' -o ' + keys[1] , ''))
            elif (len(dictSampleSeqFiles[keys][0]) == 2):
                terminal.append(bashSub("kmer_filter", dictSampleSeqFiles[keys][0], ['-1', '-2'], "--normalize " + args.depth + " --k_len " + args.k_len + ' -o ' + keys[1] , ''))

            print terminal[-1].getCommand()
            terminal[-1].runCmd("")
            sys.stderr.flush()
            #time += runSortByName.returnTime() + runView.returnTime() + htseqCmd.returnTime()
            #logFiles.append(parseOutHTseq(keys[1], keys[1].split("/")[-1]))

        #bringTogether(logFiles, os.path.join(args.finalDir, "Counts_Summary.log"))

        print "Total amount of seconds to run all samples"
        print "Seconds: " + str(time)

