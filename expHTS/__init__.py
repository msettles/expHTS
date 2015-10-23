##############################################################################
#
# expHTS is an application which organizes high-throughput sequencing
# experiments
#
##############################################################################

#from expHTS_parser import parseArgs
#from expHTS import preprocessCMD, mappingCMD, htseqCMD, parseArgs
import expHTS
from samplesheet import sampleSheet
#from preprocess import *

profile = False


def main():
    expHTS.main()

    if profile:
        #import cProfile
        #cProfile.runctx('commands[args.command].execute(args)', globals(), locals())
        return 255
    else:
        return 0;
        #commands[args.command].execute(args)

if __name__ == '__main__':
    main()
