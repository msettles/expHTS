from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOutHTseq, bringTogether
from bashSub import bashSub


#Ensures the applications are there to run
def checkPreprocessApplications():
    applications = ["samtools", "htseq-count"]
    source = ["http://samtools.sourceforge.net/", "http://www-huber.embl.de/users/anders/HTSeq/doc/install.html"]
    i = 0

    for app in applications:
        if spawn.find_executable(app) is None:
            sys.stderr.write("It doesn't look like you have app - " + app + "\n")
            sys.stderr.write("Download it here - " + source[i] + "\n")
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




class htseqCMD:

    def index(self, ref):
        if not os.path.exists(ref):
            print "Would you mind adding a gtf file? (-R) Thank you."
            exit(1)

    def execute(self, args):
        time =  0
        #Checks and makes sure user has all applications required
        checkPreprocessApplications();
        #Keeps track of log files
        logFiles = []

        # checkPreprocessApplications()
        validate = validateApp()

        #Validates and exit upon error
        validate.setValidation(True)


        #Bools for validate sample sheet
        afterPreprocess = True;
        bamFiles = True;

        dictSampleSeqFiles = validate.validateSampleSheet(args.readFolder, args.finalDir, args.samplesFile, args.force, afterPreprocess, bamFiles)

        #make sures the user added a .gtf references
        self.index(args.refGTF)



        for keys in dictSampleSeqFiles.keys():
            #checks and make sures the directories are created
            check_dir(args.finalDir)
            check_dir(keys[1])

            #Input and output files setup
            bamFile = os.path.join(keys[0], keys[0].split("/")[-1]) + ".bam"
            outFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".out"
            countFile = os.path.join(keys[1], keys[0].split("/")[-1]) + ".counts"

            # runSortByName = bashSub("samtools view -bF 0x100", [bamFile], [''], "| samtools sort -n - " + os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid"), '/dev/null')
            # htseq needs a lot of extra stuff ran on it to work. . 
            runSortByName = bashSub("samtools sort -n", [bamFile], [''],  os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid"), '/dev/null')
            print runSortByName.getCommand()
            runSortByName.runCmd("")


            #more setup type stuff for htseq
            #runSortByName = bashSub("samtools sort -n ", [os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid")], [''], os.path.join(keys[1], keys[1].split('/')[-1]) + ".byreadid" , '/dev/null')
            #print runSortByName.getCommand()
            #runSortByName.runCmd("")

            #again even more stuff for htseq
            runView = bashSub("samtools view -F 0x100 ", [os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid.bam")], [''], "> " + os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid.sam"), '/dev/null')
            print runView.getCommand()
            runView.runCmd("")


            #run htseq count
            cmdString = "htseq-count -f sam -s " + args.stranded + " -a " + args.minQual + " -t " + args.type + " -i " + args.idattr + " -m " + args.mode + " " + os.path.join(keys[1], keys[1].split('/')[-1] + ".byreadid.sam ") + args.refGTF + " 2>" + outFile + " >" + countFile


            #I'm bashsub a bit different than I should, but it still works
            htseqCmd = bashSub(cmdString, [''], [''], '', '')
            
            
            print htseqCmd.getCommand()
            htseqCmd.runCmd("")

            #Makes sure everything is pushed to the screen
            sys.stderr.flush()
            #Adds time up
            time += runSortByName.returnTime() + runView.returnTime() + htseqCmd.returnTime()

            #All count log files
            logFiles.append(parseOutHTseq(keys[1], keys[1].split("/")[-1]))
        #Creates the R summary table
        bringTogether(logFiles, os.path.join(args.finalDir, "Counts_Summary.log"))

        print "Total amount of seconds to run all samples"
        print "Seconds: " + str(time)

