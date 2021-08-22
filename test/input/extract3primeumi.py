import pysam

def extract3PrimeUmi(Infile, OutFile,  UmiTag = "MI", UmiLength = 6):
    base_complements = ''.maketrans({'A':'T', 'C':'G', 'G':'C', 'T':'A'})
    inStream = pysam.AlignmentFile(Infile, mode="rb", check_sq=False)
    OutStream = pysam.AlignmentFile(OutFile, mode="wb", template=inStream)

    for r in inStream.fetch(until_eof=True):
        xt_tag = [i for i in filter(lambda x: x[0] == 'XT', r.tags)]
        if len(xt_tag) > 0:
            # If FLAG bit 0x10 is set then SEQ is reverse complement.
            if r.is_reverse:
                s = r.seq[::-1].translate(base_complements)
                q = r.qual[::-1]
            else:
                s = r.seq
                q = r.qual
            xt_pos = xt_tag[0][1]
            umi = s[xt_pos - 1:xt_pos - 1 + UmiLength]
            umi_qual = q[xt_pos - 1:xt_pos - 1 + UmiLength]
            # RX: UMI (possibly corrected), QX: quality score for RX
            # OX: original UMI, BZ quality for original UMI
            # From SAM V1 spec:
            # The MI tag should be used for comparing UMIs. 
            # The RX tag may be used in its absence but is not guaranteed 
            # to be unique across multiple libraries.
            r.tags = r.tags + [(UmiTag, umi)]
        OutStream.write(r)

    OutStream.close()
    inStream.close()
    return

if __name__ == '__main__':
    directory = '/Users/robertshear/temp/cromwell-executions/NETseq/079825f9-8f62-414b-a36a-4c01a8ec7551/call-StarAlign/execution/'
    extract3PrimeUmi(Infile=directory+'xwt1.withXTtag.bam', OutFile=directory+'e3.bam')
