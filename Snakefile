sample_names, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")

import re                                  # Add the re import declaration to use regex
reg = re.compile('HG0...5.*')                    # Compile the regex
selected_samples = list(filter(reg.search, sample_names)) 

rule check_samples:
    input : expand("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz", sample=selected_samples)
    output: "output.txt"
    shell: "python first.py {input} {output}"
