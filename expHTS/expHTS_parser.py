# parser.py

import argparse

profile = False

# version 0.0.1 - Pre-alpha

version_num = "v0.0.1-03162015"


##############################################################################
# Preprocess app
def preprocessApp(subparsers):
    """
    preprocessApp parser parameters
    """
    #
    # Parse options
    #
    preprocess_parser = subparsers.add_parser('preprocess',
        help='Preprocess raw data high-throughput sequencing experiment')
    preprocess_parser.add_argument('-S', '--sample_metadata',
        help='file with sample metadata', action='store', type=str,
        dest='samples_file', metavar='FILENAME', default=None)
    preprocess_parser.add_argument('--debug',
        help='show traceback on error', action='store_true',
        dest="debug", default=False)
    return preprocess_parser


#####################################################################################
#  Master parser arguments
def parseArgs():
    """
    generate main parser
    """
    parser = argparse.ArgumentParser(
        description='expHTS, a python package for preprocessing of high-throughput sequencing experiments',
        epilog='For questions or comments, please contact Matt Settles <msettles@uidaho.edu>', add_help=True)
    parser.add_argument('--version', action='version', version="%(prog)s Version " + version_num)

    subparsers = parser.add_subparsers(help='commands', dest='command')

    preprocessApp(subparsers)

    args = parser.parse_args()

    return args
