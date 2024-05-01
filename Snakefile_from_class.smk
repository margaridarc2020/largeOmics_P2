
# Configuration

# Note - you can also do this through the snakemake configuration option
# read here: https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html
genome_db = "/staging/leuven/stg_00079/teaching/hg38_9/chr9.fa"
snpeff_jar = "/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/share/snpeff-5.2-0/snpEff.jar"
snpeff_genome = 'hg38'
snpeff_db_folder = '/staging/leuven/stg_00079/teaching/snpeff_db'


# list of sample names to process.
sample_names, = glob_wildcards("000.fastq/{sample}.fastq")

# ALL - the first rule is the rule snakemake automatically executes
# Note - technically you only need to add the 'tips' of the graph.
# In our case that would be the final snpeff annotated vcf file /and/
# the fastqc output - all the intermediate files will be added automatically
# however - as I develop a workflow - I tend to add files with every step that
# I add, and add them progressively here - it does not harm either.
rule all:
    input:
        fastqc_zip=expand("010.fastqc/{sample}_fastqc.zip", sample=sample_names),
        summary_png=expand("010.fastqc/{sample}_fastqc/summary.png", sample=sample_names),
        rep1=expand("010.fastqc/{sample}_fastqc/Images/per_base_quality.png", sample=sample_names),
        vcf="030.samtools/snps.vcf",
        cleaned_vcf="040.cleaned/snps.cleaned.vcf",
        snpeff = "050.snpeff/snps.annotated.vcf",



# Note - I added a number of report files here, they are automatically added to the report
#   upon rendering.
# Note 2 - although I do add the summary.txt to the output as well - it unfortunately does
#   not quite render in the final report file, but it is available as a download from within
#   the report.html. If anybody finds a solution to this, then please.
rule fastqc:
    input:
        fq="000.fastq/{file}.fastq",
    output:
        fastqc_zip="010.fastqc/{file}_fastqc.zip",
        html="010.fastqc/{file}_fastqc.html",
        summarydata="010.fastqc/{file}_fastqc/fastqc_data.txt",
        rep1=report("010.fastqc/{file}_fastqc/Images/per_base_quality.png", category="Fastqc",
                    subcategory="Per base quality", labels={"sample": "{file}"}),
        rep2=report("010.fastqc/{file}_fastqc/Images/per_base_sequence_content.png", category="Fastqc",
                     subcategory="Per base sequence content", labels={"sample": "{file}"}),
        rep3=report("010.fastqc/{file}_fastqc/summary.txt", category="Fastqc",
                    subcategory="Summary text", labels={"sample": "{file}"}),
    shell:
        """
        echo "Input Fastq: {input.fq} "
        fastqc -o 010.fastqc {input.fq} --extract

        ## Q&D testing
        ## A simple tests would be to see if output files exists
        ##   - but that is already done by snakemake

        # So - I'm testing here to see if are there FAILs in the snakemake output. If so
        # I happy to crash the workflow.

        if grep FAIL {output.rep3}; then
            # Found a fail! -
            echo "FAILED!"
            # false yields a non-zero return code - which is an error
            false
        fi

        """


# This rule is very much like the fastq rule. One input file mapped directly on two output files each.
rule bwa:
    input:
        fq="000.fastq/{sample}.fastq",
    output:
        bam = "020.bwa/{sample}.bam",
        bai = "020.bwa/{sample}.bam.bai",
    params:
        db = genome_db,
    shell:
        """
        bwa mem {genome_db} {input.fq} \
            | samtools sort - \
            > {output.bam}
        samtools index {output.bam}
        """

# This is slightly more tricky - here we have multiple BAM files, that need to
# be mapped together into one VCF file in the snp call.
# If you check the input - you see that we `expand` both generated bam files
# but have only one output file. Snakemake automatically assumes that given only
# one output file, it only needs to run this once, with both input files at the
# same time. To see exactly what it does, I've added an echo statement:
# note - the python list is printed with spaces inbetween - which is what we
# need!
rule variant_calling:
    input:
        db=genome_db,
        bams=expand("020.bwa/{sample}.bam", sample=sample_names),
    output:
        vcf="030.samtools/snps.vcf",
    shell:
        """
        echo '+-------------------------'
        echo '| processing {input.bams}'
        echo '+-------------------------'

        bcftools mpileup -Ou -f {input.db} {input.bams} \
             | bcftools call -mv -Ov -o {output.vcf}

        """


# Here we are back to 1:1 mapping of in & output
rule variant_cleanup:
    input:
        db=genome_db,
        vcf="030.samtools/snps.vcf"
    output:
        vcf="040.cleaned/snps.cleaned.vcf"
    shell:
        """
        ( cat {input.vcf} \
           | vt decompose - \
           | vt normalize -n -r {input.db} - \
           | vt uniq - \
           | vt view -f "QUAL>20" -h - \
           > {output.vcf} )


        """


rule snpeff:
    input:
        vcf = "040.cleaned/snps.cleaned.vcf",
    params:
        snpeff_db_folder = snpeff_db_folder,
        snpeff_jar = snpeff_jar,
        snpeff_genome = snpeff_genome,
    log:
        err="050.snpeff/snakemake.err",
    output:
        vcf = "050.snpeff/snps.annotated.vcf",
        html = "050.snpeff/snpEff_summary.html",
        genetxt = "050.snpeff/snpEff_genes.txt",
    shell:
        """

        mkdir -p 050.snpeff

        java -Xmx4096m -jar \
            {params.snpeff_jar} eff {params.snpeff_genome} \
            -dataDir {params.snpeff_db_folder} \
            {input.vcf} > {output.vcf}

        # move output files to the snpeff output folder
        mv snpEff_genes.txt snpEff_summary.html 050.snpeff

        """


# Here I have a python snippet that creates an image from output data
# from fastqc to include in the report
rule fastqc_report_image:
    input:
        summarytxt = "010.fastqc/{file}_fastqc/summary.txt"
    output:
        statuspng = report("010.fastqc/{file}_fastqc/summary.png",
                         category='Fastqc',
                         subcategory='Status',
                         labels={"sample": "{file}"})

    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        #load data
        data = pd.read_csv(input.summarytxt, sep="\t", header=None)
        data.columns = ['status', 'test', 'sample']

        #assign dummy x value for scatterplot
        data['x'] = 1

        #create image
        fig = plt.figure(figsize=(4,5))
        ax = plt.gca()
        sns.scatterplot(data, x='x', y='test', hue='status', s=200, ax=ax)
        ax.get_xaxis().set_visible(False)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.title(wildcards.file)
        plt.savefig(output.statuspng)
