##############################################################################
#
# expHTS is an application which organizes high-throughput sequencing
# experiments
#
##############################################################################

import sys
import os
import argparse

from expHTS_parser import parseArgs
from preprocess_app import preprocessCMD

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
