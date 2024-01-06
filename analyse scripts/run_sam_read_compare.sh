# chr1
# declare -a configs=("chr1__125000000_0_100_0_0_0_0_0_320000_0_")
# sherman_config="hg38__100000_0_100_0_1_0_0_0__0.01"

# old
# declare -a configs=("hg38__0_100000_100_0_0_0_0_0_4000000_0_" "hg38__100000_0_100_0_0_0.01_0_0_4000000_0_" "hg38__100000_0_100_0_0_0.03_0_0_4000000_0_" "hg38__100000_0_100_0_0_0.05_0_0_4000000_0_" "hg38__100000_0_100_0_0_0_0.001_0.001_
# sherman_config="hg38__100000_0_100_0_1_____0.01_"
# paired_configs=("hg38__100000_0_150_1_0_0_0_0_4000000_0_")

#configs=("chr1__10000000_0_300_0_0_0_0_0_4000000_0_")
configs=("hg38__1000000_0_100_0_0_0_0_0_4000000_0_" "hg38__1000000_0_300_0_0_0_0_0_4000000_0_" "hg38__1000000_0_100_1_0_0_0_0_4000000_0_" "hg38__1000000_0_300_1_0_0_0_0_4000000_0_" "hg38__1000000_0_100_0_1_____0.01_")
#configs=("chr21__5000000_0_100_0_0_0_0_0_40000_0_") #"chr21__5000000_0_300_0_0_0_0_0_40000_0_") #"chr21__14000001_0_100_0_0_0_0_0_400000_0_") #chr1__75000000_0_100_0_0_0.05_0.05_0.05_400000_0_") #chr1__75000000_0_100_0_0_0_0_0_400000_0_") #"chr21__14000000_0_100_0_0_0_0.05_0.05_40000_0_")
#configs=("hg38__100000_0_100_0_0_0_0.001_0.001_4000000_0_" "hg38__100000_0_100_0_0_0_0.01_0.01_4000000_0_" "hg38__100000_0_100_0_0_0_0.05_0.05_4000000_0_" "hg38__100000_0_100_0_0_0_0_0_40000000_0_" "hg38__100000_0_100_0_0_0_0_0_4000000_0_0.1" "hg38__100000_0_100_0_0_0_0_0_4000000_0_0.5" "hg38__100000_0_100_0_0_0_0_0_4000000_0_0.9" "hg38__50000_50000_100_0_0_0_0_0_4000000_0_" "hg38__0_100000_100_0_0_0_0_0_4000000_0_")
sherman_config="hg38__1000000_0_100_0_1_____0.01_"
paired_configs=("hg38__1000000_0_100_1_0_0_0_0_4000000_0_" "hg38__1000000_0_300_1_0_0_0_0_4000000_0_")

max_distance=1

declare -a algos=("aryana" "bsmap" "bwameth" "bismark" "bsbolt" "abismal")
#declare -a algos=("aryana" "aryana_em")

for config in "${configs[@]}"; do
  echo "$config"
  echo "$config" >methyl/analyze/"$config".txt
  for algo in ${algos[@]}; do
    echo "$algo"
    echo "$algo" >>methyl/analyze/"$config".txt
    if [[ $config == "$sherman_config" ]]; then
      python sam_read_compare.py ./methyl/reads/read-"$config".fastq ./methyl/outputs/result-"$config"/"$algo"/ ${max_distance} 1 >>methyl/analyze/"$config".txt
    elif [[ " ${paired_configs[*]} " =~ ${config} ]]; then
      python sam_read_compare.py -p ./methyl/reads/read-"${config}"_1.fastq ./methyl/reads/read-"${config}"_2.fastq ./methyl/outputs/result-"$config"/"$algo"/ ${max_distance} >>methyl/analyze/"$config".txt
    else
      python sam_read_compare.py ./methyl/reads/read-"$config".fastq ./methyl/outputs/result-"$config"/"$algo"/ ${max_distance} >>methyl/analyze/"$config".txt
    fi
  done
  {
    grep read methyl/analyze/"$config".txt | awk '{print $11}' | paste -s -d,
    grep read methyl/analyze/"$config".txt | awk '{print $13}' | paste -s -d,
    grep read methyl/analyze/"$config".txt | awk '{print $15}' | paste -s -d,
    grep read methyl/analyze/"$config".txt | awk '{print $17}' | paste -s -d,
  } >methyl/analyze/"$config".txt.combined
done
