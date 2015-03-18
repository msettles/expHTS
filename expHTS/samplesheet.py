
#!/usr/bin/env python

# Copyright 2015, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# samples file should have at miniumum 2 columns [SAMPLE_ID, SEQUENCE_ID]
# read in sample sheet and check for minimum required parameters
# 
#SAMPLE_ID SEQUENCE_ID
#
#LargeSeqLib   LargeSeqLib
#SmallSeqLib   SmallSeqLib

import re
import sys
import os
import os.path

class sampleSheet:

    """
    Class to read in and hold sample table information associated with Illumina sequence reads
    """
    def __init__(self, samplefile='samples.txt', sampleID='SAMPLE_ID', sequenceID='SEQUENCE_ID'):
        """
        Initialize a new sampleSheet object with the file sample list, parses and stores the sample information with
        its corresponding sequences. Class assumes the input samplefile contains the following 2 columns
        'SAMPLE_ID','SEQUENCE_ID',(defined in the header) others columns in the file are allowed and ignored
        """
        self.sampleCount = 0
        self.samples = []
        self.sequences = []
        try:
            sfile = open(samplefile, 'r')
        except IOError:
            sys.stderr.write('ERROR:[Sample sheet] cannot open\n', sfile)
            raise
        f = sfile.next()  # read first line of the file
        header = f.rstrip()
        vheader = header.split('\t')
        try:
            sampleID_index = vheader.index(sampleID)
        except ValueError:
            sys.stderr.write('ERROR:[Sample sheet] Column %s was not found in the sample sheet\n' % sampleID)
            raise
        try:
            sequenceID_index = vheader.index(sequenceID)
        except ValueError:
            sys.stderr.write('ERROR:[Sample sheet] Column %s was not found in the sample sheet\n' % sequenceID)
            raise
        try:
            for row in sfile:
                if row[0] == "#" or row[0] == "\n":  # comment or blank line
                    continue
            row = row.rstrip()
            row = row.split('\t')
            if row[sampleID_index] == '':
                raise
            if row[sequenceID_index] == '':
                raise
            self.samples.append(row[sampleID_index])
            self.sequences.append(row[sequenceID_index])
        except KeyboardInterrupt:
            sys.stderr.write('ERROR:[Sample sheet] Unexpectedly terminated')

    def getSampleCount(self):
        """
        Get the number of samples read in from the sampleFile: assumes number lines - 1
        """
        with open('samples.txt', 'r') as f:
            num_samps = sum(1 for _ in f) - 1
        return num_samps
