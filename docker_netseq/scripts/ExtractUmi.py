#!/usr/bin/env python3
import pysam
import sys

def ExtractUmi(InputBam, OutputBam, umi_length = 6, umi_tag = 'RX'):
    # input: bam file with 3' adapter code location at XT tag
    # output: RX tag added to bam file
    # TODO: parameterize umi_length
    umi_length = 6
    infile = pysam.AlignmentFile(InputBam, mode="r", check_sq=False)
    outfile = pysam.AlignmentFile(OutputBam, mode="w", template=infile)

    for r in infile.fetch(until_eof=True):
        umi = r.seq[0:umi_length]
        umi_qual = r.qual[0:umi_length]
        # RX: UMI (possibly corrected), QX: quality score for RX
        # OX: original UMI, BZ quality for original UMI
        r.tags = r.tags + [('RX', umi)]
        outfile.write(r)

    outfile.close()
    infile.close()

    return

if __name__ == '__main__':
    # TODO Validate input
    # TODO Run statistics
    # directory = '/Users/robertshear/Downloads/execution/'
    # ExtractUmi(InputBam = directory + 'xwt-1.withXTtag.bam', OutputBam = directory + 'xwt-1.trimmed.bam')
    if len(sys.argv) < 5:
        sys.exit("4 command line arguments expected")
    ExtractUmi(InputBam=sys.argv[1], OutputBam=sys.argv[2], 
        umi_length=int(sys.argv[3]), umi_tag=sys.argv[4])
