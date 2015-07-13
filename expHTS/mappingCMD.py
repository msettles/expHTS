from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOutMapping, bringTogether
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


class mappingCMD:

    def __init__(self):
        self.metaDataFolder = "MetaData"

    def createIndex(self, ref, algorithm, forceIndex):
        if os.path.exists(ref):
            if "bowtie" in algorithm:
                if not os.path.exists(ref + ".bt2") or forceIndex):
                    createTheIndex = bashSub("bowtie2-build ", [ref], [''], '', '/dev/null')
                    print createTheIndex.getCommand()
                    createTheIndex.runCmd("")
            else:
                if not os.path.exists(ref +".sa") or forceIndex:
                    createTheIndex = bashSub("bwa index ", [ref], [''], '', '/dev/null')
                    print createTheIndex.getCommand()
                    createTheIndex.runCmd("")

        else:
            print "Doesn't seem ref - " + ref + "actually exists"
            exit(1)

    def index(self, ref, algorithm, forceIndex):
        if ref == "":
            print "Would you mind adding a reference file? (-R) Thank you."
            exit(1)
        else:
            self.createIndex(ref, algorithm, forceIndex)

    def execute(self, args):
        logFiles = []  # NOT USED SO FAR
        checkPreprocessApplications()
        validate = validateApp()
        validate.setValidation(True)
        dictSampleSeqFiles = validate.validateSampleSheet(args.readFolder, args.finalDir, args.samplesFile, args.force, True)
        time = 0
        logFiles = []
        self.index(args.refFasta, args.mapping, args.forceIndex)


        for key in dictSampleSeqFiles:
            check_dir(args.finalDir)
            check_dir(key[1])
            meta = key[1]  # NOT USED SO FAR

            fileEnding = key[1].split("/")[-1]
            endString = ' 2>/dev/null | tee >(grep "^@" >' + os.path.join(key[1], fileEnding + ".header") + ') | tee >(samtools flagstat - >' + os.path.join(key[1], fileEnding + '.flagstats') + ') | samtools view -bS - | samtools sort - ' + os.path.join(key[1], fileEnding)
            SEandPE = returnReads(dictSampleSeqFiles[key])

            if SEandPE[0] != "":
                terminalString = []
            if SEandPE[1] != "":
                terminalString = []

                terminalString.append(bashSub("bwa mem", [args.threads], ['-t'], args.refFasta + " " + SEandPE[1] + " " + SEandPE[2] + " -M " + endString, "/dev/null"))
                runIndex = bashSub("samtools index ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '', '/dev/null')
                runIdxStats = bashSub("samtools idxstats ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '> ' + os.path.join(key[1], fileEnding + ".idxstats"), '/dev/null')
                #runSortByName = bashSub("samtools sort ", [os.path.join(key[1], fileEnding + ".bam")], [''], os.path.join(key[1], fileEnding + ".byreadid"), '/dev/null')

                print "___ PE COMMANDS ___"
                print terminalString[-1].getCommand()
                terminalString[-1].runCmd("")
                print runIndex.getCommand()
                runIndex.runCmd("")
                print runIdxStats.getCommand()
                runIdxStats.runCmd("")
                #if args.sortByReadID:
                #    print runSortByName.getCommand()
                #    runSortByName.runCmd("")

                sys.stderr.flush()
                time += terminalString[-1].returnTime()
                
                logFiles.append(parseOutMapping(key[1], key[1].split("/")[-1]))

            print logFiles
            bringTogether(logFiles, os.path.join(args.finalDir, "Mapping_Summary.log"))

            # logFiles.append(parseOut(key[1], key[1].split("/")[-1]))
            # bringTogether(logFiles, os.path.join(key[1].split("/")[0], "stats.log"))
            print "Total amount of seconds to run all samples"
            print "Seconds: " + str(time)

