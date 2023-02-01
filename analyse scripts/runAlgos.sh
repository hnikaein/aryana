#! /bin/bash
home_dir="/home/groups/aryanabis/"
tools_dir=${home_dir}"tools/"

conda_py3="/opt/miniconda3/bin/"
conda_py2="/opt/miniconda3/envs/py2/bin/"
snap_path=${tools_dir}"snap/"
sherman_path=${tools_dir}"Sherman/"
bs3_path=${tools_dir}"bs3/"
walt_path=${tools_dir}"walt/bin/"
abismal_path=${tools_dir}"abismal/bin/"
aryana_path=${home_dir}"aryana/"
aryana_old_path=${home_dir}"aryana_old/"
export PATH=${aryana_path}:$PATH

base_dir=${home_dir}"methyl/"
refs_dir=${base_dir}"refs/"
reads_dir=${base_dir}"reads/"
output_dir=${base_dir}"outputs/"
thread_counts=1 # XXX
cuncurrency=5 # XXX

real_read_name_arr=()
ref_name_arr=()

paired_arr=()

read_count_arr=()
cpg_read_count_arr=()
read_length_arr=()

aryana_simul_rm_arr=()
aryana_simul_ri_arr=()
aryana_simul_rd_arr=()
aryana_simul_snp_arr=()
aryana_simul_cg_arr=()

is_sherman_arr=()
sherman_er_arr=()

#config -> 1: mem-time 2:indel-snp-... 3:real
config=1 # XXX
# aligners=("aryana" "bsmap" "bwameth" "bismark" "bsbolt" "walt" "abismal") # XXX
# aligners=("bsmap" "aryana" "bwameth" "bismark" "bsbolt" "abismal")
# aligners=("aryana" "bsmap" "bwameth")
aligners=("aryana")

if [[ ${config} -eq 1 ]]
then

real_read_name_arr+=("" "" "" "" "")
ref_name_arr+=("hg38" "hg38" "hg38" "hg38" "hg38")
paired_arr+=("0") # "0" "1" "1" "0")
read_count_arr+=("100000" "100000" "100000" "100000" "100000")
cpg_read_count_arr+=("0" "0" "0" "0" "0")
read_length_arr+=("100" "300" "100" "300" "100")
aryana_simul_rm_arr+=("0" "0" "0" "0" "")
aryana_simul_ri_arr+=("0" "0" "0" "0" "")
aryana_simul_rd_arr+=("0" "0" "0" "0" "")
aryana_simul_snp_arr+=("4000000" "4000000" "4000000" "4000000" "")
aryana_simul_cg_arr+=("" "" "" "" "")
is_sherman_arr+=("0" "0" "0" "0" "1")
sherman_er_arr+=("0" "0" "0" "0" "0.01")


elif [[ ${config} -eq 2 ]]
then

real_read_name_arr+=("" "" "" "" "" "" "" "" "")
ref_name_arr+=("hg38" "hg38" "hg38" "hg38" "hg38" "hg38" "hg38" "hg38" "hg38")
paired_arr+=("0" "0" "0" "0" "0" "0" "0" "0" "0")
read_count_arr+=("100000" "100000" "100000" "100000" "50000" "0" "100000" "100000" "100000")
cpg_read_count_arr+=("0" "0" "0" "0" "50000" "100000" "0" "0" "0")
read_length_arr+=("100" "100" "100" "100" "100" "100" "100" "100" "100")
aryana_simul_rm_arr+=("0" "0" "0" "0" "0" "0" "0" "0" "0")
aryana_simul_ri_arr+=("0" "0.001" "0.05" "0.01" "0" "0" "0" "0" "0")
aryana_simul_rd_arr+=("0" "0.001" "0.05" "0.01" "0" "0" "0" "0" "0")
aryana_simul_snp_arr+=("40000000" "4000000" "4000000" "4000000" "4000000" "4000000" "4000000" "4000000" "4000000")
aryana_simul_cg_arr+=("" "" "" "" "" "" "0.1" "0.5" "0.9")
is_sherman_arr+=("0" "0" "0" "0" "0" "0" "0" "0" "0")
sherman_er_arr+=("0" "0" "0" "0" "0" "0" "0" "0" "0")


elif [[ ${config} -eq 3 ]]
then

real_read_name_arr+=("real_read_SRR3469520" "real_read_liver_small" "real_read_placenta_small")
ref_name_arr+=("hg38" "hg38" "hg38")

paired_arr+=("1") # "0" "0")

read_count_arr+=("" "" "")
cpg_read_count_arr+=("" "" "")
read_length_arr+=("125" "200" "200")

aryana_simul_rm_arr+=("" "" "")
aryana_simul_ri_arr+=("" "" "")
aryana_simul_rd_arr+=("" "" "")
aryana_simul_snp_arr+=("" "" "")
aryana_simul_cg_arr+=("" "" "")

is_sherman_arr+=("" "" "")
sherman_er_arr+=("" "" "")


elif [[ ${config} -eq 4 ]]
then

real_read_name_arr+=("" "" "" "" "")
ref_name_arr+=("chr1" "hg38" "hg38" "hg38" "hg38")
paired_arr+=("0") # "0" "1" "1" "0")
read_count_arr+=("100000" "100000" "100000" "100000" "100000")
cpg_read_count_arr+=("0" "0" "0" "0" "0")
read_length_arr+=("100" "300" "100" "300" "100")
aryana_simul_rm_arr+=("0" "0" "0" "0" "")
aryana_simul_ri_arr+=("0" "0" "0" "0" "")
aryana_simul_rd_arr+=("0" "0" "0" "0" "")
aryana_simul_snp_arr+=("4000000" "4000000" "4000000" "4000000" "")
aryana_simul_cg_arr+=("" "" "" "" "")
is_sherman_arr+=("0" "0" "0" "0" "1")
sherman_er_arr+=("0" "0" "0" "0" "0.01")


fi

total=${#paired_arr[*]}
for (( i=0; i<=$(( $total - 1 )); i++ ))
do
	ref_name=${ref_name_arr[$i]}
	real_read_name=${real_read_name_arr[$i]}
	read_count=${read_count_arr[$i]}
	cpg_read_count=${cpg_read_count_arr[$i]}

	read_length=${read_length_arr[$i]}
	is_sherman=${is_sherman_arr[$i]}

	aryana_simul_rm=${aryana_simul_rm_arr[$i]}
	aryana_simul_ri=${aryana_simul_ri_arr[$i]}
	aryana_simul_rd=${aryana_simul_rd_arr[$i]}
	aryana_simul_snp=${aryana_simul_snp_arr[$i]}
	aryana_simul_cg=${aryana_simul_cg_arr[$i]}

	sherman_er=${sherman_er_arr[$i]}

	paired=${paired_arr[$i]}

	ref_dir=${refs_dir}"ref_${ref_name}/"

	config=${ref_name}"_${real_read_name}_${read_count}_${cpg_read_count}_${read_length}_${paired}_${is_sherman}_${aryana_simul_rm}_${aryana_simul_ri}_${aryana_simul_rd}_${aryana_simul_snp}_${sherman_er}_${aryana_simul_cg}"
	echo "runing for config: ${config}"
	config_output_dir=${output_dir}"result-${config}/"
	mkdir -p ${config_output_dir}

	reads_base_name=${reads_dir}"read-${config}"
	reads=${reads_base_name}".fastq"
	reads1=${reads_base_name}"_1.fastq"
	reads2=${reads_base_name}"_2.fastq"

	sim_paired=""
	aryana_scoring=""
	if [[ ${real_read_name} != "" ]]
	then
		reads=${reads_dir}${real_read_name}".fastq"
		reads1=${reads_dir}${real_read_name}"_1.fastq"
		reads2=${reads_dir}${real_read_name}"_2.fastq"
		aryana_scoring="ar=\\\"--mp 4 --go 6 --ge 1\\\""
	elif [[ ${is_sherman} -eq 0 ]]
	then
		ratio="--mi ${ref_dir}metilation_ratio.txt"
		[[ ${aryana_simul_cg} != "" ]] && ratio="--cg ${aryana_simul_cg} --ci ${aryana_simul_cg}"

		[[ ${paired} -eq 1 ]] && sim_paired="-P"
		simul_command="read_simul -g ${ref_dir}ref.fa -n ${read_count} --ni ${cpg_read_count} -l ${read_length} ${sim_paired} -o ${reads_base_name} -m --rm ${aryana_simul_rm} --ri ${aryana_simul_ri} --rd ${aryana_simul_rd} -s ${aryana_simul_snp} -N -p -b  -i ${ref_dir}cpg_island.txt ${ratio}"
		[ ! -f ${reads} ] && [ ! -f ${reads1} ] && echo ${simul_command} && cd ${reads_dir} && ${simul_command}

		[[ ${paired} -eq 1 ]] && sed -E "s/^(@[0-9]+)_(1|2)\|.*/\1/" ${reads1} > ${reads1}.samename
		[[ ${paired} -eq 1 ]] && sed -E "s/^(@[0-9]+)_(1|2)\|.*/\1/" ${reads2} > ${reads2}.samename
	else
		export PATH=${sherman_path}:$PATH
		
		[[ ${paired} -eq 1 ]] && sim_paired="--paired_end"
		simul_command="Sherman --genome_folder ${ref_dir} -n ${read_count} -l ${read_length} ${sim_paired} --non_directional --error_rate ${sherman_er}"
		[ ! -f ${reads} ] && [ ! -f ${reads1} ] && echo ${simul_command} && cd ${reads_dir} && ${simul_command}
		
		mv simulated.fastq ${reads} 2> /dev/null
		mv simulated_1.fastq ${reads1} 2> /dev/null
		mv simulated_2.fastq ${reads2} 2> /dev/null
	fi

	for aligner_name in ${aligners[@]}
	do
		aligner_config_output_dir="${config_output_dir}${aligner_name}/"
		aligner_config_output_base_file="${aligner_config_output_dir}output"
		aligner_config_output_file="${aligner_config_output_base_file}.sam"
		aligner_ref_dir="${ref_dir}ref_${aligner_name}/"
		aligner_ref_file="${aligner_ref_dir}ref.fa"
		location=${aligner_config_output_dir}
		extra_path=""
		extra_command=""

		rm -rf ${aligner_config_output_dir} 2> /dev/null
		mkdir -p ${aligner_config_output_dir}
		cd ${aligner_config_output_dir}

		if [[ ${aligner_name} == "aryana" ]]
		then
			input_param=${reads}
			[[ ${paired} -eq 1 ]] && input_param="${reads1} ${reads2}"
			aryana_bs_params="${aligner_ref_file} ${aligner_ref_dir} ${ref_dir}cpg_island.txt ${input_param} ${aligner_config_output_base_file} ${thread_counts} ${aryana_scoring} 2>&1"
			command="aryana_bs ${aryana_bs_params}" 
			[[ ${paired} -eq 1 ]] && command="aryana_bsp ${aryana_bs_params}" 
			extra_path=${aryana_path}
		fi

		if [[ ${aligner_name} == "aryana_old" ]]
		then
			input_param=${reads}
			[[ ${paired} -eq 1 ]] && input_param="${reads1} ${reads2}"
			aryana_bs_params="${aligner_ref_file} ${aligner_ref_dir} ${ref_dir}cpg_island.txt ${input_param} ${aligner_config_output_base_file} ${thread_counts} ${aryana_scoring} 2>&1"
			command="aryana_bs ${aryana_bs_params}" 
			[[ ${paired} -eq 1 ]] && command="aryana_bsp ${aryana_bs_params}"
			extra_path=${aryana_old_path}
		fi

		if [[ ${aligner_name} == "bismark" ]]
		then
			input_param=${reads}
			maxin=""
			[[ ${paired} -eq 1 ]] && input_param="-1 ${reads1}.samename -2 ${reads2}.samename"
			[[ ${paired} -eq 1 ]] && maxin="-X 1000"
			command="bismark ${aligner_ref_dir} ${input_param} -N 1 ${maxin} --sam -o ${aligner_config_output_dir} -B output 2>&1"
		    [[ ${paired} -eq 1 ]] && extra_command="cp ${aligner_config_output_dir}output_pe.sam ${aligner_config_output_dir}output.sam 2> /dev/null;"
			extra_path=${conda_py3}
		fi

		if [[ ${aligner_name} == "bwameth" ]]
		then
			input_param=${reads}
			[[ ${paired} -eq 1 ]] && input_param="${reads1} ${reads2}"
			command="bwameth.py --reference ${aligner_ref_file} ${input_param} -t ${thread_counts} 3>${aligner_config_output_file} 2>&1 1>&3"
			extra_path=${conda_py2}
		fi

		if [[ ${aligner_name} == "bsmap" || ${aligner_name} == "bsmapz" ]]
		then
			input_param="-a ${reads}"
			[[ ${paired} -eq 1 ]] && input_param="-a ${reads1} -b ${reads2}"
			command="${aligner_name} -d ${aligner_ref_file} ${input_param} -m 0 -w 2 -R -n 1 -L ${read_length} -p ${thread_counts} -o ${aligner_config_output_file} 2>&1"
			extra_path=${conda_py3}
		fi

		if [[ ${aligner_name} == "bsbolt" ]]
		then
			input_param="-F1 ${reads}"
			[[ ${paired} -eq 1 ]] && input_param="-F1 ${reads1}.samename -F2 ${reads2}.samename"
			command="${aligner_name} Align -DB ${aligner_ref_file}.db ${input_param} -t ${thread_counts} -O ${aligner_config_output_file} 2>&1"
			extra_path=${conda_py3}
			extra_command="samtools view -h ${aligner_config_output_dir}output.sam.bam > ${aligner_config_output_dir}output.sam 2> ${aligner_config_output_dir}output.sam.err;"
		fi

		if [[ ${aligner_name} == "walt" ]]
		then
			input_param="-r ${reads}"
			[[ ${paired} -eq 1 ]] && input_param="-1 ${reads1} -2 ${reads2}"
			command="${aligner_name} -i ${aligner_ref_file}.dbindex ${input_param} -t ${thread_counts} -o ${aligner_config_output_file} 2>&1"
			extra_path=${walt_path}
		fi

		if [[ ${aligner_name} == "abismal" ]]
		then
			input_param="${reads}"
			[[ ${paired} -eq 1 ]] && input_param="${reads1} ${reads2}"
			command="${aligner_name} -i ${aligner_ref_file}.db -a -t ${thread_counts} -o ${aligner_config_output_file} ${input_param} 2>&1"
			extra_path=${abismal_path}
		fi

		if [[ ${aligner_name} == "bs_seeker" ]]
		then
			input_param="-i ${reads}"
			[[ ${paired} -eq 1 ]] && input_param="-1 ${reads1} -2 ${reads2}"
			command="bs3-align -g ${aligner_ref_file} -d ${aligner_ref_dir} ${input_param} -f sam -o ${aligner_config_output_file} 2>&1"
			rm -rf ${bs3_path}reference_genome
			ln -s ${aligner_ref_dir} ${bs3_path}reference_genome
			location=${bs3_path}
			extra_path="${conda_py2}:${bs3_path}:${snap_path}"
		fi

		timed_command="/usr/bin/time -o ${aligner_config_output_dir}time.txt -v ${command} | tee ${aligner_config_output_dir}log.txt"
		samanalyzer_command="SamAnalyzer -s -r -b -g ${aligner_ref_file} -i ${aligner_config_output_dir}output.sam -o ${aligner_config_output_dir}output.analyzed"
		tmux_command="echo \". ~/.bashrc; export PATH=${extra_path}:$PATH; echo \\\"${timed_command}\\\" >> ${aligner_config_output_dir}log.txt;cd ${location}; ${timed_command}; ${extra_command} ${samanalyzer_command}; exit; \" | bash"

		while [[ `tmux ls | wc -l` -ge ${cuncurrency} ]]
		do
			echo "waiting..."
			sleep 60 
		done
		tmux_name=${aligner_name}_${config}
		tmux_name=${tmux_name//\./}
		echo "tmux new -d -s ${tmux_name} \"${tmux_command}\""
		tmux new -d -s ${tmux_name} "${tmux_command}"
	done
done

