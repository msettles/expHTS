##############################################################################
#
# expHTS is an application which organizes high-throughput sequencing
# experiments
#
##############################################################################

#from expHTS_parser import parseArgs
from expHTS import preprocessCMD, parseArgs
from samplesheet import sampleSheet
#from preprocess import *

profile = False


def main():
    	"""
    	Entry point for expHTS application
   	"""
    	preprocess = preprocessCMD()

    	# Dictionary of subcommands, start with preprocess add more later
    	commands = {'preprocess': preprocess}

    	args = parseArgs()

    	if profile:
        	import cProfile
        	cProfile.runctx('commands[args.command].execute(args)', globals(), locals())
        	return 255
    	else:
        	commands[args.command].execute(args)

if __name__ == '__main__':
    main()
