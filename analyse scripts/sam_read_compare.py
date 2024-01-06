from sys import argv

read_files = []
is_paired = False
if argv[1] == "-p":
    is_paired = True
    read_files.append(argv[2])
    argv = argv[2:]
read_files.append(argv[1])
output_folder = argv[2]
analyze_file = output_folder + "output.analyzed"
time_file = output_folder + "time.txt"
delta = int(argv[3])
is_sherman = argv[4] if len(argv) > 4 else None
is_samename = "bsmap" in output_folder or "bsbolt" in output_folder

reads_pos = {}

for read_file in read_files:
    with open(read_file, "r") as f:
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            if line[0] != "@":
                continue
            f.readline()
            f.readline()
            f.readline()
            if is_sherman == "1":
                read_id, remain = line.split("_", 1)
            else:
                read_id, remain = line.split("|", 1)
                if is_samename:
                    read_id = read_id.split("_")[0]
            chr, remain2 = remain.split(":")
            pos = int(remain2.split("-")[0])
            read_id = int(read_id[1:])
            if read_id not in reads_pos:
                reads_pos[read_id] = []
            reads_pos[read_id].append((chr, pos))

correct = 0
total_aligned = 0
with open(analyze_file, "r") as f:
    line = f.readline()
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        parts = line.split("\t")
        if int(parts[1]) == 0:
            continue
        if is_sherman == "1":
            read_id = int(parts[0].split("_")[0])
        else:
            read_id = int(parts[0].split("|")[0])
        read_pos = reads_pos[read_id]
        total_aligned += 1
        if read_pos[0][0] != parts[2]:
            continue
        if abs(read_pos[0][1] - int(parts[3])) <= delta:
            correct += 1
            reads_pos[read_id] = reads_pos[read_id][1:]
        elif len(read_pos) > 1 and abs(read_pos[1][1] - int(parts[3])) <= delta:
            correct += 1
            reads_pos[read_id] = reads_pos[read_id][:1]
        else:
            pass
            # print(f"invalid pos {read_id} {abs(read_pos[0][1] - int(parts[3]))}")

with open(time_file, "r") as f:
    time_dict = {v[0]: v[-1] for v in [line.strip().split(":") for line in f.readlines()]}
    user_time = float(time_dict["User time (seconds)"])
    sys_time = float(time_dict["System time (seconds)"])
    mem = int(time_dict["Maximum resident set size (kbytes)"]) / 1048576.0

fdr = (total_aligned - correct) / total_aligned if total_aligned > 0 else 0
if is_paired and is_samename:
    correct /= 2.0

time = f"{user_time + sys_time:9.2f}"
mem = f"{mem:9.2f}"

print(
    f"total reads: {len(reads_pos)}\t"
    f"total aligned reads: {total_aligned}\t"
    f"correct aligned percent: {correct / len(reads_pos) * 100:4.2f}\t"
    f"fdr: {fdr:9.4f}\t"
    f"time: {time}\t"
    f"mem: {mem}")
