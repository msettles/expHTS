from validate_app import validateApp
import os
from distutils import spawn
import sys
from parse_files import parseOut, bringTogether
from bashSub import bashSub

def checkPreprocessApplications():
        applications = ["./contaminant_screen.sh", "./extract_unmapped_reads.py", "super_deduper", "sickle", "flash2"]

        for app in applications:
                if spawn.find_executable(app) is None:
                        sys.stderr.write("It doesn't look like you have app - " + app + "\n" )
                        exit(0)
                else:
                        sys.stderr.write(app + " found\n")



def returnReads(dictSampleSeqFiles):
        SE = ""
        PE1 = ""
        PE2 = ""

        #data struct 
        # { (sampleKey, seqKey) : [[SE], [SE], [PE1, PE2], [PE1, PE2]] }
        #diving into each of the sub lists in the dictionary value key
        for e in dictSampleSeqFiles:
                #if sublist only has one elment then it is SE read
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



class preprocessCMD:
        def __init__(self):
                self.metaDataFolder = "MetaData"

        def execute(self, args):
                logFiles = []
                time = 0

                checkPreprocessApplications()
                validate = validateApp()
                validate.setValidation(True)
                dictSampleSeqFiles = validate.validateSampleSheet(args.samplesDirectory, args.finalDir, args.samplesFile, args.force, False)

                for key in dictSampleSeqFiles:
                        check_dir(args.finalDir)
                        check_dir(key[1])
                        meta =  key[1]

                        SEandPE = returnReads(dictSampleSeqFiles[key])
                        #screening python scripts created in virtual enviroment
                        extract_unmapped = "python " + os.path.join(os.path.dirname(os.path.realpath(__file__)), "extract_unmapped_reads.py")
                        screen = "python " + os.path.join(os.path.dirname(os.path.realpath(__file__)), "screen.py")
                        contArgsBaseline = " -t " + args.threads
                        finalClean = "python " + os.path.join(os.path.dirname(os.path.realpath(__file__)), "cleanupWrapper.py")


                        if SEandPE[0] != "":
				terminalString = []
                                if args.contaminateFolder != "":
                                        contArgsBaseline = "-c " + args.contaminateFolder + contArgsBaseline


                                mapper = bashSub(screen, [SEandPE[0]], ['-U'], contArgsBaseline, "/dev/null")
                                cFilter = bashSub(extract_unmapped, mapper.processSub(), [''], " -o stdout" , os.path.join(meta, "SE_filter_info.log"))

				if args.skipDup == False:
                                	deduper = bashSub("super_deduper", cFilter.processSub(), ['-U'], "-p stdout", os.path.join(meta, "SE_deduper_info.log"))

                                sickleArgs = " -o " + os.path.join(key[1], "SE_not_merged.fastq")  +  " -t sanger -l " + args.minLength
                                if args.polyTrim:
                                        sickleArgs += " -a "

                                scythe = bashSub("scythe",  [args.adapter], ["-a"], deduper.processSub()[0] + " -q sanger", os.path.join(meta, "SE_scythe_info.log"))
                                sickle =  bashSub("sickle se", scythe.processSub(), ['-f'], sickleArgs, os.path.join(meta, "SE_sickle_info.log"))

                                print "___ SE COMMANDS ____"
                                print sickle.getCommand()
                                sickle.runCmd("")
                                time += sickle.returnTime()
                        if SEandPE[1] != "":
				terminalString = []
                                if args.contaminateFolder != "":
                                        contArgsBaseline = "-c " + args.contaminateFolder +  contArgsBaseline

                                terminalString.append(bashSub(screen, [SEandPE[1], SEandPE[2]], ['-1', '-2'], contArgsBaseline, "/dev/null"))
                                terminalString.append(bashSub(extract_unmapped, terminalString[-1].processSub(), [''], " -o stdout" , os.path.join(meta, "PE_filter_info.log")))
	
				if args.skipDup == False:
	                                terminalString.append(bashSub("super_deduper", terminalString[-1].processSub(), ['-i'], "-p stdout", os.path.join(meta, "PE_deduper_info.log")))


                                sickleArgs = " -m stdout -s /dev/null -t sanger -T "
                                if args.polyTrim:
                                        sickleArgs += " -a "

                                terminalString.append(bashSub("sickle pe", terminalString[-1].processSub(), ['-c'], sickleArgs , os.path.join(meta, "PE_sickle_info.log")))

                                #flash = bashSub("flash2", sickle.processSub(), ['--interleaved-input'], " -M " + args.overlapFlash + " --allow-outies -o " + key[1].split('/')[-1] + " -d " + key[1] + " 2>" + os.path.join(meta, "flash_info.log"), os.path.join(meta, "flash_info.log"))
                                #stats = bashSub("stats", flash.processSub(), [], '', os.path.join(meta, "stats.log"));
          			if args.skipFlash == False:
		                      terminalString.append(bashSub("flash2", terminalString[-1].processSub(), ['-Ti'], " -M " + args.overlapFlash + " --allow-outies -o " + key[1].split('/')[1] + " -d " + key[1] + " -To -c ", os.path.join(meta, "flash_info.log")))

                                terminalString.append(bashSub(finalClean, terminalString[-1].processSub(), [''],  " " +  str(int(args.polyTrim)) + " " + str(int(args.forceSplit)) + " " + args.minLength + " " + os.path.join(key[1], key[1].split('/')[1]), ""))



                                print "___ PE COMMANDS ___"
                                print terminalString[-1].getCommand()
                                terminalString[-1].runCmd("")
                                sys.stderr.flush()
                                time += terminalString[-1].returnTime()
                                logFiles.append(parseOut(key[1], key[1].split("/")[-1]))

                bringTogether(logFiles, os.path.join(args.finalDir, "Preprocessing_Summary.log"))
                print "Total amount of seconds to run all samples"
                print "Seconds: " + str(time)

                self.clean()


        def clean(self):
                import glob
                for f in glob.glob(".screening_cont*"):
                        os.remove(f)


