import pysam

InputFilename = '/Users/robertshear/Projects/netseq-pipeline/test/input/xwt-1.fastq.gz'
OutputFilename = '/Users/robertshear/temp/xwt-1.u.bam'

infile = pysam.AlignmentFile(InputFilename, "rb")

for s in infile:
    print(s)