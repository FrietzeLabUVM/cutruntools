# CutRunTools Installations

CutRunTools requires Python 2.7, R 3.3, Java, and Slurm job submission environment. It is designed to run on a cluster set up.

Installation also requires GCC to compile some C-based source code. 

## Prerequisites

The following tools should be already installed, and check the corresponding website to see how to install them if not. Special notes for Atactk and UCSC-tools, see below.

In the bracket is the version we have. CutRunTools may work with a lower version of each tool. 

* Trimmomatic (0.36)
* Bowtie2 (2.2.9)
* Samtools (1.3.1)
* Picard (2.8.0)
* MACS2 (2.1.1)
* MEME (4.12.0)
* Bedops (2.4.30)
* Bedtools (2.26.0)
* Atactk - we provide special install instructions since we need to patch a source file.
* UCSC-tools - we provide special install instructions.

### Atactk

This is a python2 package that determines the pA-MNase enzyme cut frequency matrix. The enzyme cuts around TF binding sites, creating DNA fragments. Frequency of enzyme cuts is thus determined by the start and end positions of DNA fragments.

The original implementation contains a mistake in computing the ending position of a fragment: it is off by 1bp. We fixed this problem. To apply the patch containing this fix and install the patched version of the package, please read `atactk.install.sh`. Make sure you have the files `make_cut_matrix.patch` and `metrics.py.patch`. Then install by:

```
source atactk.install.sh
```
This will install atactk to the user's home directory (~/.local/bin).


### UCSC-tools

We only need two commands of this package: bedgraph2BigWig and fetchChromSizes. We provide a script to automatically download and install them from UCSC genome browser:
```
source ucsc-tools.install
```
This will download the two executables in the current directory.


### Configuration file

The configuration file tells CutRunTools where to locate the prerequisite tools. This is a JSON file. It looks like below (which shows a sample `config.json`). Filling in the information should be pretty easy: in most cases we need to provide the path to the `bin` directory of each tool.

```json
{
	"Rscriptbin": "/n/app/R/3.3.3/bin",
	"pythonbin": "/n/app/python/2.7.12/bin",
	"trimmomaticbin": "/n/app/trimmomatic/0.36/bin",
	"trimmomaticjarfile": "trimmomatic-0.36.jar",
	"bowtie2bin": "/n/app/bowtie2/2.2.9/bin",
	"samtoolsbin": "/n/app/samtools/1.3.1/bin",
	"adapterpath": "/home/qz64/cutrun_pipeline/adapters", 
	"picardbin": "/n/app/picard/2.8.0/bin",
	"picardjarfile": "picard-2.8.0.jar",
	"macs2bin": "/n/app/macs2/2.1.1.20160309/bin",
	"macs2pythonlib": "/n/app/macs2/2.1.1.20160309/lib/python2.7/site-packages",
	"kseqbin": "/home/qz64/cutrun_pipeline", 
	"memebin": "/n/app/meme/4.12.0/bin", 
	"bedopsbin": "/n/app/bedops/2.4.30", 
	"bedtoolsbin": "/n/app/bedtools/2.26.0/bin",
	"makecutmatrixbin": "/home/qz64/.local/bin",
	"bt2idx": "/n/groups/shared_databases/bowtie2_indexes",
	"genome_sequence": "/home/qz64/chrom.hg19/hg19.fa",
	"extratoolsbin": "/home/qz64/cutrun_pipeline", 
	"extrasettings": "/home/qz64/cutrun_pipeline", 
	"input/output": {
		"fastq_directory": "/n/scratch2/qz64/Nan_18_aug23/Nan_run_19",
		"workdir": "/n/scratch2/qz64/workdir",
		"fastq_sequence_length": 42,
		"organism_build": "hg19"
	},
	"motif_finding": {
		"num_bp_from_summit": 150,
		"num_peaks": 5000,
		"total_peaks": 15000,
		"motif_scanning_pval": 0.0005,
		"num_motifs": 20
	},
	"cluster": {
		"email": "johndoe@gmail.com",
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
```
Please look at **the first 20 lines of this file** which concern the software installation. The rest concerns an actual analysis. 

To check if the paths are correct and if the softwares in these paths indeed work:
```
./validate.py config.json --ignore-input-output --software
```

### Download genome assemblies

CutRunTools requires the genome sequence of a specific organism build (such as hg19, hg38). We provide a script to download this automatically from NCBI. We specifically require repeat-**masked** genome sequence file. So far, the following **genome assemblies** are supported: hg19, hg38, mm10, mm9.
```
source assemblies.install
```

### Install CENTIPEDE

In R, run:
```
install.packages("CENTIPEDE", repos="http://R-Forge.R-project.org")
```

### Install CutRunTools

We first install a special trimmer we wrote `kseq`. This tool further trims the reads by 6nt to get rid of the problem of possible adapter run-through. To install:
```
source make_kseq_test.sh
```

Once completed, CutRunTools is installed.

CutRunTools does not need to be installed to any special location such as /usr/bin, or /usr/local/bin. 

It can be run directly from the directory containing the CutRunTools scripts.

The main program consists of `create_scripts.py`, `validate.py`, and a set of scripts in `aligned.aug10`, and in `macs2.narrow.aug18` which perform the important motif finding and footprinting analyses.

See USAGE.md for details. Briefly, a user writes a `config.json` configuration file for a new analysis. CutRunTools uses this to generate a set of slurm-based job-submission scripts, customized to the user's samples and environment. These job-submission scripts can be directly used by the user to start analyzing their Cut&Run samples.

