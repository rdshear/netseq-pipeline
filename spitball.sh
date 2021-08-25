time python3 <(cat <<CODE
import pysam
# input: bam file with 3' adapter code location at XT tag
# output: bam file hard clipped with RX tag added, reads without adapters removed
# TODO: parameterize umi_length
umi_length = 6
infile = pysam.AlignmentFile("xwt-1.withXTtag.bam", mode="rb", check_sq=False)
outfile = pysam.AlignmentFile("/dev/stdout", mode="w", template=infile)

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
CODE
 )  \
        | STAR --runMode alignReads \
            --genomeDir star_work \
            --runThreadN 3 \
            --readFilesIn  /dev/stdin \
            --readFilesType SAM SE \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix aligned/xwt-1. \
            --outSAMunmapped Within \
            --outSAMmultNmax 1 \
            --outSAMattributes All \
            --alignSJoverhangMin 1000 \
            --outFilterMismatchNmax 99
