# CutRunTools

This package contains the pipeline for conducting a CUT&RUN analysis.
The pipeline comprises of read trimming, alignment steps, motif finding steps, and finally the motif footprinting step. 

To get started, please see INSTALL.md about how to set up the pipeline.

Once the package is installed, please see USAGE.md to use it (but basic usage is provided below).


## Basic Usage

CutRunTools requires a JSON configuration file which specifies all that is needed to run an analysis. 
A sample configuration file is below. Note the important settings in bold.

	{
		"Rscriptbin": "/n/app/R/3.3.3/bin",
		"pythonbin": "/n/app/python/2.7.12/bin",
		"trimmomaticbin": "/n/app/trimmomatic/0.36/bin",
		"trimmomaticjarfile": "trimmomatic-0.36.jar",
		"bowtie2bin": "/n/app/bowtie2/2.2.9/bin",
		"samtoolsbin": "/n/app/samtools/1.3.1/bin",
		<b>"adapterpath"</b>: "/home/qz64",
		"picardbin": "/n/app/picard/2.8.0/bin",
		"picardjarfile": "picard-2.8.0.jar",
		"macs2bin": "/n/app/macs2/2.1.1.20160309/bin",
		"kseqbin": "/home/qz64",
		"memebin": "/home/qz64/meme/bin",
		"bedopsbin": "/home/qz64/bin",
		"bedtoolsbin": "/n/app/bedtools/2.26.0/bin",
		"makecutmatrixbin": "/home/qz64/.local/bin",
		<b>"bt2idx"</b>: "/n/groups/shared_databases/bowtie2_indexes",
		<b>"genome_sequence"</b>: "/home/qz64/chrom.hg19/hg19.fa",
		"extratoolsbin": "/home/qz64",
		"extrasettings": "/home/qz64",
		<b>"input/output"</b>: {
			<b>"fastq_directory"</b>: "/n/scratch2/qz64/Nan_18_aug23/Nan_run_19",
			<b>"workdir"</b>: "/n/scratch2/qz64/workdir",
			<b>"fastq_sequence_length"</b>: 40,
			<b>"organism_build"</b>: "hg19"
		},
		"motif_finding": {
			"num_bp_from_summit": 150,
			"num_peaks": 5000,
			"total_peaks": 15000,
			"motif_scanning_pval": 0.0005,
			"num_motifs": 20
		},
		<b>"cluster"</b>: {
			<b>"email"</b>: "bernardzhu@gmail.com",
			"step_alignment": {
				"queue": "short",
				"memory": 32000,
				"time_limit": "0-12:00"
			},
			"step_process_bam": {
				"queue": "short",
				"memory": 32000,
				"time_limit": "0-12:00"
			},
			"step_motif_find": {
				"queue": "short",
				"memory": 32000,
				"time_limit": "0-12:00"
			},
			"step_footprinting": {
				"queue": "short",
				"memory": 32000,
				"time_limit": "0-12:00"
			}
		}
	}

`fastq_directory` is the directory containing paired-end CUT&RUN sequences (with _R1_001.fastq.gz and _R2_001.fastq.gz suffix). `organism_build` is one of supported genome assemblies: hg38, hg19, mm10, and mm9. `adapterpath` contains Illumina Truseq3-PE adapter sequences (we provide them). `genome_sequence` is the whole-genome **masked** sequence which matches with the appropriate organism build.

### Create job submission scripts

	./create_scripts.py config.json

This creates a set of Slurm job-submission scripts based on the configuration file above (named config.json) in the `workdir` directory. The scripts can be directly, easily executed.

	./validate.py config.json

This script checks that your configuration file is correct and all paths are correct.

### Four-step process

With the scripts created, we can next perform the analysis.

1. **Read trimming, alignment.** We suppose the `workdir` is defined as `/n/scratch2/qz64/workdir`

	cd /n/scratch2/qz64/workdir
	sbatch ./integrated.sh CR_BCL11A_W9_r1_S17_R1_001.fastq.gz

The parameter is the fastq file. Always use the _R1_001 version of the pair.

2. **BAM processing, peak calling.** It marks duplicates in bam files, and filter fragments by size.

	cd aligned.aug10
	sbatch ./integrated.step2.sh CR_BCL11A_W9_r1_S17_aligned_reads.bam


3. **Motif finding.** CutRunTools uses MEME-chip for de novo motif finding on sequences surrounding the peak summits.

	cd ../macs2.narrow.aug18
	sbatch ./integrate.motif.find.sh CR_BCL11A_W9_r1_S17_aligned_reads_peaks.narrowPeak

By default, CutRunTools keeps duplicate fragments. If instead users wish to use deduplicate version, 

	cd ../<b>macs2.narrow.aug18.dedup</b>
	sbatch ./integrate.motif.find.sh CR_BCL11A_W9_r1_S17_aligned_reads_peaks.narrowPeak


4. **Motif footprinting.**

	cd ../macs2.narrow.aug18
	sbatch ./integrate.footprinting.sh CR_BCL11A_W9_r1_S17_aligned_reads_peaks.narrowPeak

Beautiful footprinting figures will be located in the directory `fimo.result`. Footprinting figures are created for every motif found by MEME-chip, but only the right motif (associated with TF) will have have a proper looking shape. Users can scan through individual motif's footprint.


