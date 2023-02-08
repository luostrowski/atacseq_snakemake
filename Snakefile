
# Snakefile to analyze ATACseq PE data

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

localrules: collect_bam_metrics

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
	    expand("data/trimming/{sample}.paired_{dir}.fq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
	    #expand("data/bowtie2/{sample}.bam", sample = SAMPLES),
        expand("data/bowtie2/a_sorted_bams/{sample}.sorted.bam", sample = SAMPLES),
        expand("data/bowtie2/a_sorted_bams/{sample}.read_count.txt", sample = SAMPLES),
        expand("data/bowtie2/a_sorted_bams/{sample}.sorted.bam.bai", sample = SAMPLES),
        expand("data/bowtie2/b_noMT/{sample}.sorted.noMT.bam", sample = SAMPLES),
        expand("data/bowtie2/b_noMT/{sample}.read_count.txt", sample = SAMPLES),
        expand("data/bowtie2/c_rmDups/{sample}.sorted.noMT.rmDups.bam", sample = SAMPLES),
        expand("data/bowtie2/c_rmDups/{sample}.dup_metrics.txt", sample = SAMPLES),
        expand("data/bowtie2/c_rmDups/{sample}.read_count.txt", sample = SAMPLES),
        expand("data/bowtie2/d_mapq30/{sample}.sorted.noMT.rmDups.mapq30.bam", sample = SAMPLES),
        expand("data/bowtie2/d_mapq30/{sample}.read_count.txt", sample = SAMPLES),
        expand("data/bowtie2/e_proper_pair/{sample}.sorted.noMT.rmDups.mapq30.proper.bam", sample = SAMPLES),
        expand("data/bowtie2/e_proper_pair/{sample}.read_count.txt", sample = SAMPLES),
        "data/bowtie2/metrics/summary_metrics.txt"

rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip",
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"


rule trimmomatic:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
	    rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/trimming/{sample}.paired_R1.fq.gz",
	    rev = "data/trimming/{sample}.paired_R2.fq.gz",
	    fwd_unpaired = "data/trimming/{sample}.unpaired_R1.fq.gz",
	    rev_unpaired = "data/trimming/{sample}.unpaired_R2.fq.gz"
    conda:
        "envs/trimmomatic.yaml"
    params:
        trimmer = ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"],
	    adapters = config["adapters"]
    shell:
        "trimmomatic PE -phred33 {input.fwd} {input.rev} {output.fwd} {output.fwd_unpaired} {output.rev} {output.rev_unpaired} {params.adapters} {params.trimmer}"

rule bowtie2:
    input:
        fwd = "data/trimming/{sample}.paired_R1.fq.gz",
	    rev = "data/trimming/{sample}.paired_R2.fq.gz"
    output:
        bam_file = "data/bowtie2/a_sorted_bams/{sample}.bam"
    conda:
        "envs/bowtie2.yaml"
    params:
        genome = config["bt2_genome"]
    shell:
        "bowtie2 -p 8 -X 2000 -x {params.genome} -1 {input.fwd} -2 {input.rev} | samtools view -bS - > {output.bam_file}"

rule samtools_sort:
    input:
        bam_file = "data/bowtie2/a_sorted_bams/{sample}.bam"
    output:
        outfile = "data/bowtie2/a_sorted_bams/{sample}.sorted.bam",
        count = "data/bowtie2/a_sorted_bams/{sample}.read_count.txt"
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        samtools sort -o {output.outfile} {input.bam_file}
        rm {input.bam_file}
        samtools view -F 0x4 {output.outfile} | cut -f 1 | sort | uniq | wc -l > {output.count}
        """

rule samtools_index:
    input:
        bam_file = "data/bowtie2/a_sorted_bams/{sample}.sorted.bam"
    output:
        bai = "data/bowtie2/a_sorted_bams/{sample}.sorted.bam.bai"
    conda:
        "envs/bowtie2.yaml"
    shell:
        "samtools index {input.bam_file} {output.bai}"

rule remove_MT:
    input:
        bam_file = "data/bowtie2/a_sorted_bams/{sample}.sorted.bam"
    output:
        outfile = "data/bowtie2/b_noMT/{sample}.sorted.noMT.bam",
        count = "data/bowtie2/b_noMT/{sample}.read_count.txt"
    conda:
        "envs/bowtie2.yaml"
    params:
        mt_chr = config["mt_chr"]
    shell:
        """
        samtools idxstats {input.bam_file} | cut -f 1 | grep -v {params.mt_chr} | xargs samtools view -b {input.bam_file} > {output.outfile}
        samtools view -F 0x4 {output.outfile} | cut -f 1 | sort | uniq | wc -l > {output.count}
        """

rule picard_mark_dups:
    input:
        bam_file = "data/bowtie2/b_noMT/{sample}.sorted.noMT.bam"
    output:
        outfile = "data/bowtie2/c_rmDups/{sample}.sorted.noMT.rmDups.bam",
        metrics = "data/bowtie2/c_rmDups/{sample}.dup_metrics.txt",
        count = "data/bowtie2/c_rmDups/{sample}.read_count.txt"
    conda:
        "envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input.bam_file} \
            O={output.outfile} \
            M={output.metrics} \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
        samtools view -F 0x4 {output.outfile} | cut -f 1 | sort | uniq | wc -l > {output.count}
        """

rule filter_mapq30:
    input:
        bam_file = "data/bowtie2/c_rmDups/{sample}.sorted.noMT.rmDups.bam"
    output:
        outfile = "data/bowtie2/d_mapq30/{sample}.sorted.noMT.rmDups.mapq30.bam",
        count = "data/bowtie2/d_mapq30/{sample}.read_count.txt"
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        samtools view -h -q 30 {input.bam_file} > {output.outfile}
        samtools view -F 0x4 {output.outfile} | cut -f 1 | sort | uniq | wc -l > {output.count}
        """

rule proper_pairs:
    input:
        bam_file = "data/bowtie2/d_mapq30/{sample}.sorted.noMT.rmDups.mapq30.bam"
    output:
        outfile = "data/bowtie2/e_proper_pair/{sample}.sorted.noMT.rmDups.mapq30.proper.bam",
        count = "data/bowtie2/e_proper_pair/{sample}.read_count.txt"
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        samtools view -b -f 0x02 {input.bam_file} > {output.outfile}
        samtools view -F 0x4 {output.outfile} | cut -f 1 | sort | uniq | wc -l > {output.count}
        """

rule collect_bam_metrics:
    input:
        expand("data/bowtie2/e_proper_pair/{sample}.read_count.txt", sample = SAMPLES),
        expand("data/bowtie2/c_rmDups/{sample}.dup_metrics.txt", sample = SAMPLES)
    output:
        #bam_metrics = "data/bowtie2/metrics/bam_metrics.txt"
        #dup_metrics = "data/bowtie2/metrics/dup_metrics.txt"
        outfile = "data/bowtie2/metrics/summary_metrics.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        bt2_outdir = "data/bowtie2/",
        bam_metrics = "data/bowtie2/metrics/bam_metrics.txt",
        dup_dir = "data/bowtie2/c_rmDups/",
        dup_metrics = "data/bowtie2/metrics/dup_metrics.txt"
    shell:
        """
        python scripts/summarize_bam_counts.py -d {params.bt2_outdir} -o {params.bam_metrics}
        python scripts/parse_picard_dups.py -d {params.dup_dir} -o {params.dup_metrics}
        join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,2.2,1.5,1.6 <(head -n 2 {params.bam_metrics} && tail -n +3 {params.bam_metrics} | sort) <(head -n 2 {params.dup_metrics} && tail -n +3 {params.dup_metrics} | sort) > {output.outfile}
        """
