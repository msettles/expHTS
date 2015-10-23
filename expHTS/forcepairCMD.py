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


class forcepairCMD:

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
        dictSampleSeqFiles = validate.validateSampleSheet(args.readFolder, args.forcepairFolder, args.samplesFile, args.force, True)

        for keys in dictSampleSeqFiles.keys():
            check_dir(args.forcepairFolder)
            check_dir(keys[1])
            terminal = []
            #countFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".counts"
            print keys
            if (len(dictSampleSeqFiles[keys][0]) == 3):
                cmdString = "cp " + dictSampleSeqFiles[keys][0][0] + " " + os.path.join(keys[1], ".") + " & "
                cmdString += " cp " + dictSampleSeqFiles[keys][0][1] + " " + os.path.join(keys[1], ".") + ";"
                r1File = os.path.join(keys[1], dictSampleSeqFiles[keys][0][0].split("/")[-1])
                r2File = os.path.join(keys[1], dictSampleSeqFiles[keys][0][1].split("/")[-1])

                cmdString += """awk '{ printf "%s/1\\n", $0
                        getline
                        print substr($0, 0, length/2)
                        getline
                        print $0
                        getline
                        print substr($0, 0, length/2)
                      }' """ + dictSampleSeqFiles[keys][0][2] + " >> " + r1File +  " & " 
                cmdString += """awk 'BEGIN {
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
                   }'  """ + dictSampleSeqFiles[keys][0][2] + " >> " + r2File

                terminal.append(bashSub(cmdString, [''], [''], '', ''))
            
                print terminal[-1].getCommand()
                terminal[-1].runCmd("")
                sys.stderr.flush()
        time += terminal[-1].returnTime()

        #logFiles.append(parseOutHTseq(keys[1], keys[1].split("/")[-1]))
        #bringTogether(logFiles, os.path.join(args.finalDir, "Counts_Summary.log"))

        print "Total amount of seconds to run all samples"
        print "Seconds: " + str(time)

