# preprocess_app.py
import sys


class preprocessCMD:
    """
    preprocessing application for expHTS experiments
    """
    def __init__(self):
        pass

    def execute(self, args):
        # ----------------------- options input files -----------------------
        if args.samples_file is None:
            sFile = None
            sys.stderr.write("No sample file identified\n")
        else:
            sFile = args.samples_file
        print sFile

        return 0
