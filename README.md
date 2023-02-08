# atacseq_snakemake
Snakemake pipeline for ATACseq PE data

To use for a new project:

  - clone this repository from github
  - symlink the raw fastq files into the data/fastq directory (ln -s /path/to/file1.fq.gz .)
  - ensure fastq file names are proper format (mv file1.fq.gz sample_R1.fastq.gz)
  - edit the proj_config.yaml file if necessary:
    - adapters: default is Nextera PE, usually works with ATACseq
    - bt2_genome: path to reference genome bowtie2 indexes    
    - mt_chr: mitochondrial genome name from reference genome, usually "chrM" or "MT"    

Example execution of snakefile:

Note, use the -n parameter first to test that the structure works and not actually execute anythting.   

snakemake --use-conda --jobs 100 --cluster-config cluster.json --cluster "sbatch --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr}"
