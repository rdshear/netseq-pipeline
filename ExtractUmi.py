import pysam
# input: bam file with 3' adapter+bar code location at XT tag
# output: bam file hard clipped with RX tag added, reads without adapters removed
umi_length = 6
base_complements = ''.maketrans({'A':'T', 'C':'G', 'G':'C', 'T':'A'})
#infile = pysam.AlignmentFile("/dev/stdin", mode="rb", check_sq=False)
infile = pysam.AlignmentFile("/Users/robertshear/temp/cromwell-executions/NETseq/fe04dfc0-6fd6-46b0-80b6-a076e9c45fee/call-StarAlign/execution/foo.bam", mode="rb", check_sq=False)
outfile = pysam.AlignmentFile("/dev/stdout", mode="wb", template=infile)

for r in infile.fetch(until_eof=True):
    xt_tag = [i for i in filter(lambda x: x[0] == 'XT', r.tags)]
    if len(xt_tag) > 0:
        s = r.seq
        q = r.qual
        xt_pos = xt_tag[0][1]
        umi = s[xt_pos - 1:xt_pos - 1 + umi_length]
        umi_qual = q[xt_pos - 1:xt_pos - 1 + umi_length]
        # RX: UMI (possibly corrected), QX: quality score for RX
        # OX: original UMI, BZ quality for original UMI
        r.seq = s[0:xt_pos]
        r.qual = q[0:xt_pos]
        r.tags = r.tags + [('RX', umi)]
        outfile.write(r)

outfile.close()
infile.close()