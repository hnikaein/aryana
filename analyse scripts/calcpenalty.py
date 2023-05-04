# ReadName        Aligned AlnChr  AlnPos  RefLength       AlnMismatch     AlnGapOpen      AlnGapExt       RealChr
# RealPos CorrectPos      RealMismatch    RealGapOpen     RealGapExt      SamSeq  AlnCIGAR        AlnSeq  RealCIGAR
# RealSeq

algos = ["aryana_old", "bsmap"]
algos_len = len(algos)
samples = [
    # "/home/groups/aryanabis/methyl/outputs/result-hg38_real_read_liver_small2.fastq___200________/",
    # "/home/groups/aryanabis/methyl/outputs/result-hg38_real_read_liver_small.fastq___200________/",
    # "/home/groups/aryanabis/methyl/outputs/result-hg38_real_read_placenta_small.fastq___200________/",
    "/home/groups/aryanabis/methyl/outputs/result-hg38_real_read_SRR19154020___150_1_______/",
]
paired_end = True
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
            if read_name not in reads:
                reads[read_name] = [[], []]  # * algos_len
            if int(parts[1]) == 0:
                continue
            elif int(parts[1]) != 1:
                print("!!!!!!!!! -> ", parts)
            mis = int(parts[5])
            go = int(parts[6])
            ge = int(parts[7])
            mat = int(parts[4]) - mis - go - ge
            pen = 4 * mis + 6 * go + 1 * ge - 4 * mat
            if len(reads[read_name][algo_i]) > 0 and reads[read_name][algo_i][0][1] != parts[2]:
                continue
            reads[read_name][algo_i].append((pen, parts[2]))
        f.close()
    unique_alignment = [0, 0]  # * algos_len
    for read_name, penalties in reads.items():
        if paired_end:
            for i in range(2):
                if len(penalties[i]) < 2:
                    penalties[i] = []
        diff_len = len(penalties[0]) - len(penalties[1])
        if diff_len > 0:
            unique_alignment[0] += 2
        elif diff_len < 0:
            unique_alignment[1] += 2
        else:
            penalties[0].sort()
            penalties[1].sort()
            for i in range(min(len(penalties[0]), len(penalties[1]))):
                dif = penalties[1][i][0] - penalties[0][i][0]
                if dif < -900:
                    print(dif, read_name)
                if dif not in differention:
                    differention[dif] = 1
                else:
                    differention[dif] += 1
    print(sorted(differention.items()))
    print(unique_alignment)
