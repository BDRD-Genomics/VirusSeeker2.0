#!/bin/bash
shopt -s extglob
basedir=/path/to/VS/dir
execdir="$basedir/scripts/"
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#source /path/to/miniforge3/etc/profile.d/conda.sh
#conda activate vs

timestamp() {
       	date +%Y-%m-%d_%H%M%S
}

display_usage() { 
	echo -e "Usage:\n$0 [-d host_db for host removal] [-h help] [-a assembly_type] [-m memory in GB] [-o output_directory] [-t threads] [-v run in Virome mode] <-r readlist file>\n" >&2
}

while getopts "h?d:m:o:r:t:va:k:" opt; do
       	case "$opt" in
		d)	host_db=$OPTARG ;; # host removal (path to host fasta file)
	       	h|\?)	display_usage 
			exit 1
			;;
	       	m)  	slurm_mem=$OPTARG ;; # memory per process in GB
	       	o)  	outdir=$basedir/VS_output/$OPTARG ;; # e.g. batch1
	       	r)  	readlist=$OPTARG ;; # file which contains path-to-reads
	       	t)  	slurm_cpus_per_task=$OPTARG ;; # threads per stage e.g. 64
	       	a)		assembly=$OPTARG ;; #assembly type (s for illumina PE, l for ONT long reads, h for hybrid assembly)
	       	v)  	mode="V" ;; # run in Virome mode, defaults to Discovery mode
	       	k)      assembly_mode=$OPTARG ;; #assembly type (m for metagenomics, i for isolate)
       	esac
done

if ((OPTIND == 1))
then 
	echo "No options specified" >&2
	display_usage
	exit 1
fi


if [[ -z ${readlist} ]]
then
	echo "-r readlist required" >&2
	exit 1
fi

[[ -z ${host_db} ]] && host_db=0
[[ -z ${mode} ]] && mode="D"
[[ -z ${assembly} ]] && assembly="s"
[[ -z ${assembly_mode} ]] && assembly_mode="m"
[[ -z ${slurm_mem} ]] && slurm_mem=0 # defaults to value in VS.config
[[ -z ${slurm_cpus_per_task} ]] && slurm_cpus_per_task=0 # defaults to value in VS.config
[[ -z ${outdir} ]] && outdir=$basedir/VS_output/out.$(timestamp) && mkdir ${outdir}
[ ! -d ${outdir} ] && mkdir -p ${outdir}
log="$outdir/queue.log"
exec >>${log} 2>&1

echo -e "\n*** Starting virusseeker queue process at $(timestamp) ***"
#echo -e "${assembly}"
count_inputs=()

for sample in $(cat ${readlist})
do
	sname=$(basename ${sample} | \
		sed -n 's:_S[0-9]*_L00[0-9]_R[1-2]_001.fastq.gz::p;s:_S[0-9]*_R[1-2]_001.fastq.gz::p;s:_R[1-2]_001.fastq.gz::p;s:_[Rr][1-2].fastq.gz::p;s:_ONT.fastq.gz::p;s:_ont.fastq.gz::p;s:.fastq.gz::p' )
	[ ! -d ${outdir}/${sname} ] && mkdir -p ${outdir}/${sname}
	count_inputs+="${sname} "
	
	if [[ "$assembly" == "s" ]]
	then
		r1=(${outdir}/${sname}/${sname}_*[rR]1*)
		r2=(${outdir}/${sname}/${sname}_*[rR]2*)
		if [[ ! -s ${r1} ]] || [[ ! -s ${r2} ]]
			then
			echo -e "Copying ${sample} at $(timestamp)\n"
			cp -rv ${sample} ${outdir}/${sname}
		fi

		#if [[ ${sample} = *[rR]2* && "${assembly}" == "s" ]]
		if [ $(grep -ow ${sname} <<< ${count_inputs[A]} | wc -l) -eq 2 ]
		# 
		then
			touch ${outdir}/${sname}/${sname}_adapter.txt
			ln -sf ${outdir}/${sname}/${sname}_*R1* \
					${outdir}/${sname}/${sname}_SE1.fastq.gz 
			ln -sf ${outdir}/${sname}/${sname}_*R2* \
						${outdir}/${sname}/${sname}_SE2.fastq.gz 
			echo -e "Queueing ${sname} at $(timestamp)\n"
			perl ${execdir}/runVS.pl ${outdir}/${sname} ${host_db} ${slurm_cpus_per_task} ${slurm_mem} ${assembly} ${assembly_mode} 0  0 ${mode}; ;; #removes human 
		fi

	elif [[ "${assembly}" == "l" ]]
	then
		r=(${outdir}/${sname}/${sname}.fastq.gz)
		if [[ ! -s ${r} ]]
			then
			echo -e "Copying ${sample} at $(timestamp)\n"
			cp -rv ${sample} ${outdir}/${sname}/${sname}.fastq.gz
		fi

		if [[ -s ${r} ]]
		then
			#only one read so doesn't matter which sample you're on
			echo -e "wrong path"
			touch ${outdir}/${sname}/${sname}_adapter.txt
			ln -sf ${outdir}/${sname}/${sname} \
						${outdir}/${sname}/${sname}_LR.fastq.gz 
			echo -e "Queueing ${sname} at $(timestamp)\n"
			perl ${execdir}/runVS.pl ${outdir}/${sname} ${host_db} ${slurm_cpus_per_task} ${slurm_mem} ${assembly} ${assembly_mode} 0  0 ${mode}; ;; #removes human 
		fi
	elif [[ "$assembly" == "h" ]]
	then
		r=(${outdir}/${sname}/${sname}.fastq.gz)
		r1=(${outdir}/${sname}/${sname}_*[rR]1*)
		r2=(${outdir}/${sname}/${sname}_*[rR]2*)

		if [[ ! -s ${r1}  || ! -s ${r2}  || ! -s ${r} ]]
		then
			echo -e "Copying ${sample} at $(timestamp)"
			cp -rv ${sample} ${outdir}/${sname}
		fi
		
		echo ""
		echo "$count_inputs"
		echo -e $(grep -o $sname <<< ${count_inputs[*]} | wc -l)
		if [ $(grep -ow $sname <<< ${count_inputs[*]} | wc -l) -eq 3 ]
		then
			#echo -e "All there!"
			touch ${outdir}/${sname}/${sname}_adapter.txt
			ln -sf ${outdir}/${sname}/${sname}_*R1* \
					${outdir}/${sname}/${sname}_SE1.fastq.gz 
			ln -sf ${outdir}/${sname}/${sname}_*R2* \
					${outdir}/${sname}/${sname}_SE2.fastq.gz 
			ln -sf ${outdir}/${sname}/${sname}.fastq.gz \
					${outdir}/${sname}/${sname}_LR.fastq.gz 
			echo -e "Queueing ${sname} at $(timestamp)\n"
			perl ${execdir}/runVS.pl ${outdir}/${sname} ${host_db} ${slurm_cpus_per_task} ${slurm_mem} ${assembly} ${assembly_mode} 0  0 ${mode}; ;; #removes human 
		fi
	fi
done
echo -e "Input array = ${count_inputs[@]}"
echo -e "*** Finished at $(timestamp) ***"
exit 0
