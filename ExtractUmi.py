import pysam

umi_length = 6
base_complements = ''.maketrans({'A':'T', 'C':'G', 'G':'C', 'T':'A'})
infile = pysam.AlignmentFile("/dev/stdin", mode="rb")
outfile = pysam.AlignmentFile("/dev/stdout", mode="wb", template=infile)

for r in infile.fetch(until_eof=True):
    xt_tag = [i for i in filter(lambda x: x[0] == 'XT', r.tags)]
    if len(xt_tag) > 0:
        # STAR translates seq to reverse complement, restore it to get true UMI
        if r.is_reverse:
            s = r.seq[::-1].translate(base_complements)
            q = r.qual[::-1]
        else:
            s = r.seq
            q = r.qual
        xt_pos = xt_tag[0][1]
        umi = s[xt_pos - 1:xt_pos - 1 + umi_length]
        umi_qual = q[xt_pos - 1:xt_pos - 1 + umi_length]
        # RX: UMI (possibly corrected), QX: quality score for RX
        # OX: original UMI, BZ quality for original UMI
        r.tags = r.tags + [('RX', umi)]
        # Only generate records with UMI's
        outfile.write(r)

outfile.close()
infile.close()
