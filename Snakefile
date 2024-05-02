# run in the terminal: export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/bin/:$PATH

sample_names, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")
genome_db = "/lustre1/project/stg_00079/teaching/hg38_21/chr21.fa"

import re                                  # Add the re import declaration to use regex
reg = re.compile('HG0...5.*')                    # Compile the regex
selected_samples = list(filter(reg.search, sample_names)) 

rule all:
    input: 
        "output.txt", 
        expand("qc_output_folder/fastqc_files.fastqc/{sample}_fastqc.html", sample=selected_samples),
        expand("bwa_files.bwa/{sample}.bam", sample=selected_samples),
        expand("bwa_files.bwa/{sample}.bam.bai", sample=selected_samples),
        "vcf_files.samtools/snps.vcf"

rule fastqc:
    input:
        fastq_files="/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz",
    output:
        report="qc_output_folder/fastqc_files.fastqc/{sample}_fastqc.html"
    shell:
        """
        mkdir -p qc_output_folder
        fastqc -o qc_output_folder/fastqc_files.fastqc {input.fastq_files} --extract
        ls -l qc_output_folder/fastqc_files.fastqc
        echo {output.report}
        ls {output.report}
        """

rule check_samples:
    input : expand("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz", sample=selected_samples)
    output: "output.txt"
    shell: "python first.py {input} {output}"


rule bwa:
    input:
        fastq_files ="/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz",
    output:
        bam = "bwa_files.bwa/{sample}.bam",
        bai = "bwa_files.bwa/{sample}.bam.bai",
    params:
        db = genome_db,
    shell:
        """
        bwa mem {params.db} {input.fastq_files} \
            | samtools sort - \
            > {output.bam}
        samtools index {output.bam}
        """

rule variant_calling:
    input:
        db=genome_db,
        bams=expand("bwa_files.bwa/{sample}.bam", sample=selected_samples),
    output:
        vcf="vcf_files.samtools/snps.vcf",
    shell:
        """
        bcftools mpileup -Ou -f {input.db} {input.bams} \
             | bcftools call -mv -Ov -o {output.vcf}
        """

