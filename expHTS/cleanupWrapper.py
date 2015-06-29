from finalCleanup import finalCleanup
import sys
import os
def main():
	print '/'.join(sys.argv[5].split('/')[:-1])
	finalCleanup(int(sys.argv[2]),   int(sys.argv[3]), int(sys.argv[4]), os.path.join('/'.join(sys.argv[5].split('/')[:-1]), "logtestFile.log"),  str(sys.argv[1]), str(sys.argv[5]) + "_R1.fastq", str(sys.argv[5]) + "_R2.fastq",  str(sys.argv[5]) + "_SE.fastq")


main() 
