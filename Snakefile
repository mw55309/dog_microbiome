
def get_accessions():
	
	accns = list()
	dogfile = open("dogs.txt", mode="r")

	for row in dogfile:
		accns.append(row.strip('\n').split('\t')[1])

	return(list(set(accns)))


def get_dogs():
	
	dogs = list()
	dogfile = open("dogs.txt", mode='r')

	for row in dogfile:
		dogs.append(row.strip('\n').split('\t')[1])

	return(list(set(dogs)))


def get_r1_accessions_for_dog(wildcards):

	dogs = list()
	dogfile = open("dogs.txt", mode="r")
	
	mydog = wildcards.id

	for row in dogfile:
		arr = row.strip('\n').split('\t')
		
		if arr[1] == mydog:
			dogs.append("trimmed/" + arr[0] + "_1" + ".t.fastq.gz")

	dogs = sorted(list(set(dogs)))
	return(dogs)


def get_r2_accessions_for_dog(wildcards):

	dogs = list()
	dogfile = open("dogs.txt", mode="r")

	mydog = wildcards.id

	for row in dogfile:
		arr = row.strip('\n').split('\t')

		if arr[1] == mydog:
			dogs.append("trimmed/" + arr[0] + "_2" + ".t.fastq.gz")

	dogs = sorted(list(set(dogs)))
	return(dogs)
	


DOGS = get_dogs()
ACCN = get_accessions()

print(DOGS)

######################################################
#
# 
# rule "all" is the default rule that Snakemake runs
# this rule basically pulls data through the entire
# pipeline by specifying the final outputs of the
# pipeline as input. The rule does nothing
#
#
######################################################

rule all:
	input: expand("megahit/{sample}", sample=DOGS)


rule megahit:
	input:
		R1="combined/{id}.R1.fq.gz",
		R2="combined/{id}.R2.fq.gz"
	output:
		di=directory("megahit/{id}")
	conda: "envs/megahit.yaml"
	threads: 8
	resources:
		mem_mb=64000, disk_mb=80000
	shell: 
		'''
		mkdir -p {output.di} && megahit --continue --k-list 31,59,87 --kmin-1pass -m 0.95 --min-contig-len 1500 -m 64000000000 -t {threads} -1 {input.R1} -2 {input.R2} -o {output.di} && rm -rf {output.di}/intermediate_contigs
		'''


rule combine:
	input: 
		R1=get_r1_accessions_for_dog,
		R2=get_r2_accessions_for_dog
	output: 
		R1="combined/{id}.R1.fq.gz",
		R2="combined/{id}.R2.fq.gz"
	resources:
		mem_mb=8000, disk_mb=50000
	shell:
		'''
		cat {input.R1} > {output.R1} && cat {input.R2} > {output.R2} && rm {input.R1} {input.R2}
		'''


rule download_and_trim:
	output:
		R1="trimmed/{id}_1.t.fastq.gz",
		R2="trimmed/{id}_2.t.fastq.gz"
	params:
		id="{id}"
	conda: "envs/cutadapt.yaml"
	threads: 2
	resources:
		mem_mb=16000, disk_mb=18000
	shell: 
		'''
		fastq-dump --split-files --split-spot --gzip {params.id} && cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o {output.R1} -p {output.R2} -O 5 --minimum-length=50 {params.id}_1.fastq.gz {params.id}_2.fastq.gz
		'''
