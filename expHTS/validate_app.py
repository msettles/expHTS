import sys
import os

class validateApp:
    
    sampleFiles = {}

    def __init__(self):
        self.exitValidate = False
        self.debug = True
        self.verbose = True


    def validateSampleSheet(self, dirSample, finalDir, sampleSheet, force, afterPreprocess, bamFiles):

        linenum = 0
        
        #Test to see if sample directory exists
        if not os.path.exists(dirSample):
            self.exitTime("Directory " + dirSample + " is not found")
        #See if there is a sample sheet
        if os.path.exists(sampleSheet):
            f = open(sampleSheet, "r");
        #Is there just nothing
        else:
            self.exitTime("No sample sheet named " + sampleSheet + " was found")


        #reads sample sheet
        #there are 3 booleans that could be a bit confusing
        #force will overwrite the current info in the directories
        #preprocess is a special case - previous to preprocess the sampleID and seqID are different
        #therefore we need to read in column 1 and column 2, after proprocess, we just need col 2
        #htseq takes in back files, not fastq, so that has a special function
        for e in f.readlines():
            if linenum != 0:
                self.sampleSequenceID(dirSample, finalDir, e.split("\t"), force, afterPreprocess, bamFiles)

            linenum += 1


        
        for key in self.sampleFiles:
            self.sampleFiles[key].sort(key=lambda x: len(x))

        return self.sampleFiles


    #Gathers data form sample sheet
    def sampleSequenceID(self, dirSample, finalDir, seqID, force, afterPreprocess, bamFiles):
        
        #comment, ignore line
        if seqID[0][0] == "#":
            pass
        elif len(seqID) >= 2:

            #as stated earlier, preprocess needs both columns
            if (afterPreprocess == False):
                seqID[0] = os.path.join(dirSample, seqID[0].rstrip()) 
                seqID[1] = os.path.join(finalDir, seqID[1].rstrip())
                #Test samples
                self.finalDirTest(seqID[1], force)
                self.directoryFiles(seqID)
            #after preprocess, only the second column is needed
            else:
                seqID[0] = os.path.join(dirSample, seqID[1].rstrip()) 
                seqID[1] = os.path.join(finalDir, seqID[1].rstrip())
                self.finalDirTest(seqID[1], force)

                if bamFiles == False:
                    #If fastq files are being looked for
                    self.directoryFiles(seqID)
                else:
                    #Some applications will be looking for a bam file (such as htseq)
                    self.directoryFilesBam(seqID, seqID[1].split('/')[-1]);    
        else:
            self.exitTime("There wasn't two columns in the sample file file")



    
    def finalDirTest(self, sampleID, force):
        #If there was data all ready there - just keep it and don't touch it unless force is set
        if os.path.exists(sampleID) and not force:
            self.exitTime(sampleID + " was all ready created. Use the -w or --overwrite option to overwrite")
        elif os.path.exists(sampleID):
            #Just shoot a nice warning letting the user know there is no hope
            #and that there data is now being overwritten
            print "Warning"
            print "Overwrite was turned on - overwriting " + sampleID + "\n"
            
            

    #sets up the directory dictionary
    #set up key with tuple (sample and seq)
    def directoryFilesBam(self, sampleSeq, fileName):
        sampleSeq = tuple(sampleSeq)
        #print sampleSeq
        directoryTest = sampleSeq[0].rstrip();
        #fastqCount insures at least one fastq files is under the directory
        bamCount = 0


        #Going down the directory
        if self.testDirectory(directoryTest):
            for subdir, dir, files in os.walk(directoryTest):
                for file in files:
                    file = os.path.abspath(os.path.join(directoryTest, file))
                    #checks for bam format
                    if ".bam" in file and fileName in file and not ".bai" in file:
                        bamCount += 1
                        #see if the array has been initialized yet
                        if not sampleSeq in self.sampleFiles:
                            self.sampleFiles[sampleSeq] = []
                        
                        #adds to the dictionary (no file should be named the same so this is okay)
                        self.sampleFiles[sampleSeq].append(file)

            if (bamCount == 0):
                self.exitTime("No bam files were found under this directory - " + directoryTest)

        else:
            self.exitTime("Directory " + directoryTest + " does not exists")
            



    #THIS WILL NEED TO BE CHANGED
    #Currently, we are just using the Casava(?) Format that dictates *_R1.fastq and *_R2.fastq
    #A future enhancement where each file is checked and queued to ensure Pairs, SE, and Interleaved
    #are all account for.
    def directoryFiles(self, sampleSeq):
        sampleSeq = tuple(sampleSeq)
        #print sampleSeq
        directoryTest = sampleSeq[0].rstrip();
        #fastqCount insures at least one fastq files is under the directory
        fastqCount = 0
        if self.testDirectory(directoryTest):
            for subdir, dir, files in os.walk(directoryTest):
                for file in files:
                    file = os.path.abspath(os.path.join(directoryTest, file))

                    #Yes, this is only checking for fastq(.gz) files that have _R1 in them
                    #Previous to preprocess, _R1 alone will be a SE read
                    if "_R1" in file and ".fastq" in file:
                        fastqCount += 1

                        #starts an empty array with the sample sheet
                        if not sampleSeq in self.sampleFiles:
                            self.sampleFiles[sampleSeq] = []
                        
                        self.sampleFiles[sampleSeq].append(self.isPairedEnd(file))


            if (fastqCount == 0):
                self.exitTime("No fastq files were found under this directory - " + directoryTest)

        else:
            self.exitTime("Directory " + directoryTest + " does not exists")
            
    def exitTime(self, exitString):
        print exitString
        if self.exitValidate:
            print "Failure in Validation - exiting process"
            exit(-1)
        else:
            print "Failure in Validation - going to continue to continue to validate sample sheet"

    #returns 2 elements in list if paired end reads or 1 element in list if it is a SE
    def isPairedEnd(self, fileRead1):

        file2 = fileRead1.replace("_R1", "_R2")
        SE = fileRead1.replace("_R1", "_SE")

        if os.path.exists(SE) and os.path.exists(file2):
            return [fileRead1, file2, SE]
        if os.path.exists(file2):
            return [fileRead1, file2]
        else:
            return [fileRead1]
        

    #test if directory exists . . . meh. probably can get ride of it
    def testDirectory(self, directoryTest):
        return os.path.isdir(directoryTest)

    #outputs the single end paired end reads    
    def infoOutput(self):
        #print self.sampleFiles
        for key in self.sampleFiles:
            #print "SEQUENCE " + key[0] + " SAMPLE " + key[1]
            for files in self.sampleFiles[key]:
                if len(files) == 1:
                    pass
                    #print "SE Files: " + str(files).strip("[]")
                else:
                    pass
                    #print "PE Files: " + str(files).strip("[]")

    def setValidation(self, exitOut):
        self.exitValidate = exitOut

        
    #tuple key (seq ID - sample ID)
    #value - 1 value is SE reads and 2 values is PE
    def dictionaryFilesReturn(self):
        return self.sampleFiles


