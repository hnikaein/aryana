# ReadName        Aligned AlnChr  AlnPos  RefLength       AlnMismatch     AlnGapOpen      AlnGapExt       RealChr
# RealPos CorrectPos      RealMismatch    RealGapOpen     RealGapExt      SamSeq  AlnCIGAR        AlnSeq  RealCIGAR
# RealSeq

algos = ["aryana", "bwameth"]
algos_len = len(algos)
samples = [
    #    "/home/groups/aryanabis/methyl/outputs/result-hg38_real_read_liver_small2.fastq___200________/",
    "/home/groups/aryanabis/methyl/outputs/result-hg38_real_read_liver_small.fastq___200________/",
    "/home/groups/aryanabis/methyl/outputs/result-hg38_real_read_placenta_small.fastq___200________/",
]
for sample in samples:
    reads = {}
    differention = {}
    for algo_i, algo in enumerate(algos):
        fname = sample + algo + "/output.analyzed"
        print(fname)
        f = open(fname, "r")
        f.readline()
        while True:
            parts = f.readline().strip().split("\t")
            if len(parts) < 7:
                break
            read_name = parts[0]
            if not read_name in reads:
                reads[read_name] = [-1] * algos_len
            if int(parts[1]) == 0:
                continue
            elif int(parts[1]) != 1:
                print("!!!!!!!!! -> ", parts)
            mis = int(parts[5])
            go = int(parts[6])
            ge = int(parts[7])
            mat = int(parts[4]) - mis - go - ge
            pen = 4 * mis + 6 * go + 1 * ge - 4 * mat
            reads[read_name][algo_i] = pen
        f.close()
    unique_alignment = [0] * algos_len
    for read_name, penalties in reads.items():
        for i in range(algos_len):
            if penalties[i] == -1:
                for j in range(algos_len):
                    if penalties[j] != -1:
                        unique_alignment[j] += 1
                break
        else:
            dif = penalties[1] - penalties[0]
            if dif not in differention:
                differention[dif] = 1
            else:
                differention[dif] += 1
    print(sorted(differention.items()))
    print(unique_alignment)
