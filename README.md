# ATAC-seq snakemake pipeline

Snakemake pipeline for processing ATACseq PE data.

## How to use
### 1. Prepare work environment

- Clone this repository in the working directory

```
git clone https://github.com/luostrowski/atacseq_snakemake.git
```

- Navigate into the data/fastq directory

```
cd data/fastq
```

- Symlink files to the directory

```
ln -s /absolute/path/to/file_R1.fq.gz .
ln -s /absolute/path/to/file_R2.fq.gz .
```

- Rename symlink to proper format

```
mv file_R1.fq.gz sample_R1.fastq.gz
mv file_R2.fq.gz sample_R2.fastq.gz
```

### 2. Adjust configuration file

The `proj_config.yaml` file contains the wildcards that the script will use. 

- Change the path to the reference genome (it must be indexed).

- Change the path to the annotation gtf file. 

- Adjust treat samples if necessary.

- Change trimmomatic adapters if necessary. Default is Nextera PE, usually works with ATAC-seq.

- Change the mitochondrial genome name (`mt_chr`) if necessary. Options are: "chrM" or "MT". 

- Adjust `effective_genome_size` for deeptools bamcoverage:   
      - GRCh37: 2864785220    
      - GRCh38: 2913022398    
      - GRCm37: 2620345972    
      - GRCm38: 2652783500    

- Adjust other parameters if necessary:
    - `macs2_genome`: two letter code for MACS2 to recognize genome   
    - `contrasts`: desired contrasts for differential analysis   
    - `analysis_method`: for diffbind, either DBA_DESEQ2 or DBA_EDGER   
    - `min_overlap`: for diffbind merging, min number of samples peak found to retain it   
    - `fdr_cutoff`: padj cutoff for defining significant diff peaks   
    - `data_source`: for chipseeker, genome used
    - `organism`: for chipseeker, organism latin name with underscore (no spaces)

### 4. Execute the snakefile

- First try a dry run to make sure the structure works.

```
snakemake --dry-run
```

- Adjust the `cluster.json` file to run on a cluster submitting slurm jobs.

    - The `--latency-wait 60` parameter is required when using slurm to make sure outfiles are complete.

```
snakemake --use-conda --jobs 100 --latency-wait 60 --cluster-config cluster.json --cluster "sbatch --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr}"
```

- If it is necessary to re-run the script, the incomplete files will need to be re-generated. 

```
snakemake --rerun-incomplete --use-conda --jobs 100 --latency-wait 60 --cluster-config cluster.json --cluster "sbatch --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr}"
```

## Directed acyclic graph of the steps

![DAG](dag_simplified.pdf)
