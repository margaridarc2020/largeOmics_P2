sample_names, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")
#print(sample_names)

rule first:
    input : expand("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz", sample=sample_names)
    output: "coucou.txt"
    shell: "python first.py {input} {output}"
