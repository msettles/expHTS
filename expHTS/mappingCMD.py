from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOutMapping, bringTogether
from bashSub import bashSub



def checkPreprocessApplications():
    applications = ["bwa", "samtools", "bowtie2"]
    source = ["http://bio-bwa.sourceforge.net/", "http://samtools.sourceforge.net/", "http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/"]
    i = 0

    for app in applications:
        if spawn.find_executable(app) is None:
            sys.stderr.write("It doesn't look like you have app - " + app + "\n")
            sys.stderr.write("Download it here - " + source[i] + "\n")
            exit(0)
        else:
            sys.stderr.write(app + " found\n")
        i+= 1

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
                if not os.path.exists(ref + ".bt2") or forceIndex:
                    createTheIndex = bashSub("bowtie2-build ", [ref], [''], '', '/dev/null')
                    print createTheIndex.getCommand()
                    createTheIndex.runCmd("")
            else:
                if not os.path.exists(ref + ".sa") or forceIndex:
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
        time = 0
        checkPreprocessApplications();
        logFiles = []

        # checkPreprocessApplications()
        validate = validateApp()
        validate.setValidation(True)
        dictSampleSeqFiles = validate.validateSampleSheet(args.readFolder, args.finalDir, args.samplesFile, args.force, True)

        self.index(args.refFasta, args.mapping, args.forceIndex)

        for key in dictSampleSeqFiles:
            check_dir(args.finalDir)
            check_dir(key[1])
            meta = key[1]  # NOT USED SO FAR

            fileEnding = key[1].split("/")[-1]
            endString = ' 2>/dev/null | tee >(samtools flagstat - >' + os.path.join(key[1], fileEnding + '.flagstats') + ') | samtools view -bS - | samtools sort - ' + os.path.join(key[1], fileEnding)
            SEandPE = returnReads(dictSampleSeqFiles[key])
            files = dictSampleSeqFiles[key][0]
            RGstring = "-R '@RG\tID:" + fileEnding + "\tSM:" + fileEnding + "\tPL:ILLUMINA\tLB:whatever\tPU:whatever\tDS:Paired'"

            if SEandPE[0] != "":
                terminalString = []
            if len(files) == 3 and args.forcePairs:
                awkR1 = """
                awk '{
                   printf "%s/1\\n", $0
                      getline
                         print substr($0, 0, length/2)
                              getline
                                print $0 
                                  getline
                                    print substr($0, 0, length/2)
                                     }' """ + files[2]
                awkR2 = """awk 'BEGIN {
                  j = n = split("A C G T", t)
                    for (i = 0; ++i <= n;)
                        map[t[i]] = t[j--]
                          }
                            {
                               printf "%s/2\\n", $0
                                  getline
                                     for (i = length; i > length/2 ; i--) {
                                         printf "%s", map[substr($0, i, 1)]
                                            }
                                              printf "\\n"
                                                getline
                                                  print $0 
                                                    getline
                                                      for (i = length; i > length/2; i--) {  
                                                          printf "%s", substr($0, i, 1)
                                                            }
                                                              printf "\\n"
                                                               }'  """ + files[2]
                terminalString = []
                terminalString.append(bashSub("bwa mem -M " + RGstring, [str(int(args.threads))], ['-t'], args.refFasta + " <(cat " + files[0] + " <(" +  awkR1 + ")) <(cat " + files[1] + " <(" + awkR2 + ")) "  + endString, "/dev/null"))
                runIndex=bashSub("samtools index ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '', '/dev/null')
                runIdxStats=bashSub("samtools idxstats ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '> ' + os.path.join(key[1], fileEnding + ".idxstats"), '/dev/null')


            elif len(files) == 3:
                terminalString = []
                terminalString.append(bashSub("bwa mem -M " + RGstring, [str(int(args.threads)/2)], ['-t'], args.refFasta + " " + files[0] + " " + files[1] , "/dev/null"))
                terminalString.append(bashSub("bwa mem -M " + RGstring, [str(int(args.threads)/2)], ['-t'], args.refFasta + " " + files[2] , "/dev/null"))
                terminalString.append(bashSub("samtools merge - " + terminalString[-1].processSub()[0] + " " + terminalString[-2].processSub()[0] + " " + endString , [""], [""], "", "/dev/null"))

                runIndex=bashSub("samtools index ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '', '/dev/null')
                runIdxStats=bashSub("samtools idxstats ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '> ' + os.path.join(key[1], fileEnding + ".idxstats"), '/dev/null')

            elif SEandPE[1] != "":
                RGstring = "-R '@RG\tID:" + fileEnding + "\tSM:" + fileEnding + "\tPL:ILLUMINA\tLB:whatever\tPU:whatever\tDS:Paired'"
                terminalString = []
                terminalString.append(bashSub("bwa mem -M " + RGstring, [args.threads], ['-t'], args.refFasta + " " + SEandPE[1] + " " + SEandPE[2] + endString, "/dev/null"))
                runIndex = bashSub("samtools index ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '', '/dev/null')
                runIdxStats = bashSub("samtools idxstats ",  [os.path.join(key[1], fileEnding + ".bam")], [''], '> ' + os.path.join(key[1], fileEnding + ".idxstats"), '/dev/null')


            print "___ PE COMMANDS ___"
            print terminalString[-1].getCommand()
            terminalString[-1].runCmd("")
            print runIndex.getCommand()
            runIndex.runCmd("")
            print runIdxStats.getCommand()
            runIdxStats.runCmd("")

            sys.stderr.flush()
            time += runIndex.returnTime() + runIdxStats.returnTime() + terminalString[-1].returnTime()

            logFiles.append(parseOutMapping(key[1], key[1].split("/")[-1]))

        bringTogether(logFiles, os.path.join(args.finalDir, "Mapping_Summary.log"))

        print "Total amount of seconds to run all samples"
        print "Seconds: " + str(time)

