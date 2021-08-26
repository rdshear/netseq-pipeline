#!/usr/bin/env python3
import pysam
import sys

def ExtractUmi(InputBam, OutputBam, umi_length = 6, umi_tag = 'RX', adapter_trim_tag = 'XT'):
    # input: bam file with 3' adapter code location at XT tag
    # output: bam file hard clipped with RX tag added, reads without adapters removed
    # TODO: parameterize umi_length
    umi_length = 6
    infile = pysam.AlignmentFile(InputBam, mode="rb", check_sq=False)
    outfile = pysam.AlignmentFile(OutputBam, mode="wb", template=infile)

    for r in infile.fetch(until_eof=True):
        xt_tag = [i for i in filter(lambda x: x[0] == 'XT', r.tags)]
        if len(xt_tag) > 0:
            s = r.seq
            q = r.qual
            xt_pos = xt_tag[0][1]
            umi = s[0:umi_length]
            umi_qual = q[0:umi_length]
            # TODO? Remove XT tag bcause the reads are trimmed
            # RX: UMI (possibly corrected), QX: quality score for RX
            # OX: original UMI, BZ quality for original UMI
            r.seq = s[umi_length:xt_pos]
            r.qual = q[umi_length:xt_pos]
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
    if len(sys.argv) < 6:
        sys.exit("5 command line arguments expected")
    ExtractUmi(InputBam=sys.argv[1], OutputBam=sys.argv[2], 
        umi_length=int(sys.argv[3]), umi_tag=sys.argv[4], adapter_trim_tag=sys.argv[5])