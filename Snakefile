sample_names, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")

import re                                  # Add the re import declaration to use regex
reg = re.compile('HG0...5.*')                    # Compile the regex
selected_samples = list(filter(reg.search, sample_names)) 

rule all:
    input: 
        "output.txt", 
        expand("output_folder/fastqc_files.fastqc/{sample}_fastqc.html", sample=selected_samples)

rule fastqc:
    input:
        fastq_files="/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz",
    output:
        report="output_folder/fastqc_files.fastqc/{sample}_fastqc.html"
    shell:
        """
        mkdir -p output_folder
        fastqc -o output_folder/fastqc_files.fastqc {input.fastq_files} --extract
        ls -l output_folder/fastqc_files.fastqc
        echo {output.report}
        ls {output.report}
        """

rule check_samples:
    input : expand("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz", sample=selected_samples)
    output: "output.txt"
    shell: "python first.py {input} {output}"



