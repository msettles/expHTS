#!/usr/bin/env python

'''
Extract reads which aren't mapped from a SAM or SAM.gz file.
Behavior for PE:
  -Write out PE only if both do not map (if either of the pair maps, neither is retained)
Behavior for SE:
  -Write out SE if they don't map

Iterate over a SAM or SAM.gz file. take everything where the 3rd and
4th flag bit are set to 1 and write reads out to files.

0x1 template having multiple segments in sequencing
0x2 each segment properly aligned according to the aligner
0x4 segment unmapped
0x8 next segment in the template unmapped
0x10 SEQ being reverse complemented
0x20 SEQ of the next segment in the template being reversed
0x40 the first segment in the template
0x80 the last segment in the template
0x100 secondary alignment
0x200 not passing quality controls
0x400 PCR or optical duplicate

TODO:
1) Add support for retaining both reads if one of a pair don't map but the other does
2) Add support for retaining the pair (or SE) if a read maps with low mapq

Note:
  It is necessary to double check that both pairs of the PE read really exist in the SAM
file just in case it somehow gets disordered. This is taken care of by keeping the PE
reads in a set of dictionaries and then deleting them once the pair is written.
In the case where a read is somehow labeled as paired, but the pair doesn't exist, the
read is NOT written.
'''
import sys
import os
from optparse import OptionParser  # http://docs.python.org/library/optparse.html
import gzip


usage = "usage: %prog [options] -o output_base inputfile.SAM"
parser = OptionParser(usage=usage, version="%prog 2.0.1")
parser.add_option('-u', '--uncompressed', help="leave output files uncompressed",
                  action="store_true", dest="uncompressed")
parser.add_option('-o', '--output_base', help="output file basename",
                  action="store", type="str", dest="output_base", default="screened")
parser.add_option('-t', '--tab-seperated', help="seperated out the output in tab",
                  action="store_true",  dest="tabSeperated")
parser.add_option('-v', '--verbose', help="verbose output",
                  action="store_false", dest="verbose", default=True)

(options,  args) = parser.parse_args()  # uncomment this line for command line support

if len(args) == 1:
    infile = args[0]
    # Start opening input/output files:
    if not os.path.exists(infile):
        print >> sys.stderr, "Error, can't find input file %s" % infile
        sys.exit()

    if infile.split(".")[-1] == "gz":
        insam = gzip.open(infile, 'rb')
    else:
        insam = open(infile, 'r')
else:
    # reading from stdin
    insam = sys.stdin

base = options.output_base

PE1 = {}
PE2 = {}
contig_map = {}
interleaved = False


def writeread(ID, r1, r2):
    if interleaved:
        if options.tabSeperated is True:
            # read1
            print ID + "\t" + r1[0] + "\t" + r1[1] + "\t" + r2[0] + "\t" + r2[1] + "\n"
        else:
            print "@" + ID + "#0/1"
            print r1[0]
            print '+\n' + r1[1]
            # read2
            print "@" + ID + "#0/2"
            print r2[0]
            print '+\n' + r2[1]
    else:
        # read1
        outPE1.write("@" + ID + "#0/1" '\n')
        outPE1.write(r1[0] + '\n')
        outPE1.write('+\n' + r1[1] + '\n')
        # read2
        outPE2.write("@" + ID + "#0/2" '\n')
        outPE2.write(r2[0] + '\n')
        outPE2.write('+\n' + r2[1] + '\n')

i = 0
PE_written = 0
SE_written = 0
SE_open = False
PE_open = False
line2 = []

for line in insam:
    # Comment/header lines start with @
    if line[0] != "@" and len(line.strip().split()) > 2:
        line2 = line.strip().split()
        flag = int(line2[1])

        if (flag & 0x100):  # secondary alignment
            continue

        i += 1
        # Handle SE:
        # unapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
        if (flag & 0x1 == 0) and (flag & 0x4):
            ID = line2[0].split("#")[0]
            if not SE_open:
                if base == "stdout":
                    interleaved = True
                elif options.uncompressed:
                    outSE = open(base + "_SE.fastq", 'w')
                else:
                    outSE = gzip.open(base + "_SE.fastq.gz", 'wb')
                SE_open = True
            # interleaved just means to stdout in this case
            if (interleaved):
                if options.tabSeperated is True:
                    print ID + "\t" + line2[9] + "\t" + line2[10] + "\n"
                else:
                    print "@" + ID
                    print line2[9]
                    print '+\n' + line2[10]
            else:
                outSE.write("@" + ID + '\n')
                outSE.write(line2[9] + '\n')
                outSE.write('+\n' + line2[10] + '\n')
            SE_written += 1
            continue

        # Handle PE:
        # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped, 0x80 the last segment in the template
        if ((flag & 0x1) and (flag & 0x4) and (flag & 0x8)):
            if not PE_open:
                if base == "stdout":
                    interleaved = True
                elif options.uncompressed:
                    outPE1 = open(base + "_PE1.fastq", 'w')
                    outPE2 = open(base + "_PE2.fastq", 'w')
                else:
                    outPE1 = gzip.open(base + "_PE1.fastq.gz", 'wb')
                    outPE2 = gzip.open(base + "_PE2.fastq.gz", 'wb')
                PE_open = True
            if (flag & 0x40):  # is this PE1 (first segment in template)
                # PE1 read, check that PE2 is in dict and write out
                ID = line2[0].split("#")[0]
                r1 = [line2[9], line2[10]]  # sequence + qual
                if ID in PE2:
                    writeread(ID, r1, PE2[ID])
                    del PE2[ID]
                    PE_written += 1
                else:
                    PE1[ID] = r1
                continue
            elif (flag & 0x80):  # is this PE2 (last segment in template)
                # PE2 read, check that PE1 is in dict and write out
                ID = line2[0].split("#")[0]
                r2 = [line2[9], line2[10]]
                if ID in PE1:
                    writeread(ID, PE1[ID], r2)
                    del PE1[ID]
                    PE_written += 1
                else:
                    PE2[ID] = r2
                continue
        # was mapped, count it up
        # if line2 != []:
        #     contig = line2[2]
        #     if contig in contig_map.keys():
        #         if (flag & 0x1 == 0):  # SE
        #             contig_map[contig]["SE"] += 1
        #         elif (flag & 0x40):  # PE, Just count the first in the pair
        #             contig_map[contig]["PE"] += 1
        #     else:
        #         contig_map[contig] = {}
        #         if (flag & 0x1 == 0):  # SE
        #             contig_map[contig]["SE"] = 1
        #             contig_map[contig]["PE"] = 0
        #         elif (flag & 0x40):  # PE, Just count the first in the pair
        #             contig_map[contig]["SE"] = 0
        #             contig_map[contig]["PE"] = 1

# for k in contig_map.keys():
#     print >> sys.stderr,  "\tFound %s: percent: %.2f, PE mapped: %s, SE mapped: %s" % (k, (2*PE_written+SE_written)/i, contig_map[k]["PE"], contig_map[k]["SE"])

print >> sys.stderr, "Records processed: %s | PE_written: %s | SE_written: %s | Discarded: %s " % (i, PE_written, SE_written, i-(PE_written*2+SE_written))

if base != "stdout":
    if PE_open:
        outPE1.close()
        outPE2.close()
    if SE_open:
        outSE.close()
