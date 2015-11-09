#!/usr/bin/python


from subprocess import PIPE, STDOUT, Popen, check_call
import os
import time
import sys
import signal
import optparse
import time

#Class I use to do bash substitution easier
#To initialize, (execute, [files], [file args], additional args, output_file (/dev/null maybe)
class bashSub:
    app = ""
    files = ""
    file_args = ""    
    args = ""
    f_info = ""    
    cmd = ""

    def __init__(self, app, files, file_args, args, info_file):
        self.app = app
        self.files = files
        self.file_args = file_args
        self.args = args
        self.f_info = info_file
        self.time = 0

        #Insure file and file args are the same
        if len(files) != len(file_args):
            print >> sys.stderr, "Error: File args need to be the same length"

        self.cmd = app
        
        #Turns array into string
        for i in range(0, len(files)):
            self.cmd += " " + file_args[i] + " " + files[i] + " "
        
        self.cmd += " " + args

        #sets to output file            
        if self.f_info != "":
            self.cmd += " 2>" + self.f_info
        
    #just returns command string
    def getCommand(self):
        return self.cmd



    def runCmd(self, additional_args):
        start = time.time()
        p2 = Popen(self.cmd + additional_args, stdout=sys.stdout, stdin=PIPE, stderr=sys.stderr, shell=True, executable = "/bin/bash");
        (out, error) = p2.communicate()
        sys.stderr.flush();

        if p2.poll() != 0:
            print "-----------------Command Failed---------------"
            print self.cmd + additional_args
            print "-----------------Error message----------------"
            newcmd = (self.cmd + additional_args).split(">")[0][:-1]
            p2 = Popen(newcmd, stdout=PIPE, stdin=PIPE, stderr=PIPE, shell=True, executable = "/bin/bash");
            (out, error) = p2.communicate()
            print error
            sys.exit();
            
        end = time.time()
        self.time = end - start

    def returnTime(self):
        return self.time


    def processSub(self):
        return ["<(" + self.cmd + ")"]

