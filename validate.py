#!/usr/bin/python
import shutil
import sys
import os
import re
import json

def check_program_exists(path, program):
	if not os.path.isfile(path + "/" + program):
		print program, "is not found in", path
		return False
	return True
	

if __name__=="__main__":
	f = open("config.json")
	config = json.load(f)
	f.close()

	#validate the input files
	flist = set([])
	for filename in os.listdir(config["input/output"]["fastq_directory"]):
		if filename.endswith("fastq.gz"):
			flist.add(filename)

	for x in flist:
		m1 = re.match("(.*)_R1_001.fastq.gz", x)
		m2 = re.match("(.*)_R2_001.fastq.gz", x)
		if m1 is None and m2 is None:
			print "Error:", x, " no pattern of _R1_001.fastq.gz, or _R2_001.fastq.gz detected"
			sys.exit(1)
		if m1 is not None:
			x2 = m1.group(1) + "_R2_001.fastq.gz"
			if not x2 in flist:
				print x, "does not have a corresponding _R2_001.fastq.gz file"
				sys.exit(1)
		if m2 is not None:
			x2 = m2.group(1) + "_R1_001.fastq.gz"
			if not x2 in flist:
				print x, "does not have a corresponding _R1_001.fastq.gz file"
				sys.exit(1)

	if not check_program_exists(config["bowtie2bin"], "bowtie2"):
		sys.exit(1)
	if not check_program_exists(config["memebin"], "meme-chip"):
		sys.exit(1)

	if not check_program_exists(config["bedopsbin"], "bedops"):
		sys.exit(1)
	if not check_program_exists(config["bedopsbin"], "gff2bed"):
		sys.exit(1)
	if not check_program_exists(config["bedopsbin"], "sort-bed"):
		sys.exit(1)

	if not check_program_exists(config["samtoolsbin"], "samtools"):
		sys.exit(1)
	if not check_program_exists(config["trimmomaticbin"], config["trimmomaticjarfile"]):
		sys.exit(1)
	if not check_program_exists(config["macs2bin"], "macs2"):
		sys.exit(1)
	if not check_program_exists(config["picardbin"], config["picardjarfile"]):
		sys.exit(1)
	if not check_program_exists(config["kseqbin"], "kseq_test"):
		sys.exit(1)
	if not check_program_exists(config["bedtoolsbin"], "bedtools"):
		sys.exit(1)
	if not check_program_exists(config["makecutmatrixbin"], "make_cut_matrix"):
		sys.exit(1)
	if not check_program_exists(config["Rscriptbin"], "Rscript"):
		sys.exit(1)
	if not check_program_exists(config["pythonbin"], "python"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "bedGraphToBigWig"):
		sys.exit(1)
	if not check_program_exists(config["extratoolsbin"], "fetchChromSizes"):
		sys.exit(1)

	
	if not check_program_exists(config["extrasettings"], "filter_below.awk"):
		sys.exit(1)

	if not check_program_exists(config["adapterpath"], "Truseq3.PE.fa"):
		sys.exit(1)

	if not check_program_exists(config["extrasettings"], "%s.chrom.sizes" % 
	config["input/output"]["organism_build"]):
		sys.exit(1)

	org = config["input/output"]["organism_build"]
	if org=="hg38":
		org = "GRCh38"
	for ff in ["%s.1.bt2" % org, "%s.2.bt2" % org, "%s.3.bt2" % org, \
	"%s.4.bt2" % org, "%s.rev.1.bt2" % org, "%s.rev.2.bt2" % org]:	
		if not check_program_exists(config["bt2idx"], ff):
			sys.exit(1)
			
