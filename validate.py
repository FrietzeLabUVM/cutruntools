#!/usr/bin/python
import shutil
import sys
import os
import re
import json
import argparse

def check_program_exists(path, program):
	if not os.path.isfile(path + "/" + program):
		print program, "is not found in", path
		return False
	return True
	

if __name__=="__main__":
	parser = argparse.ArgumentParser(prog="validate")
	parser.add_argument("--software", action="store_true")
	parser.add_argument("--ignore-input-output", action="store_true")
	parser.add_argument("config", type=file)

	args = vars(parser.parse_args(sys.argv[1:]))
	#f = open(args["config"]) #config.json
	f = args["config"]
	config = json.load(f)
	f.close()

	#check input and output portion of config
	if not args["ignore_input_output"]:
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

	#valide software path
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

	#check genome sequence exists, and check bowtie2 indices
	if not args["ignore_input_output"]:
		if not check_program_exists(os.path.dirname(config["genome_sequence"]), "%s.chrom.sizes" % 
		config["input/output"]["organism_build"]):
			sys.exit(1)

		org = config["input/output"]["organism_build"]
		if org=="hg38":
			org = "GRCh38"
		for ff in ["%s.1.bt2" % org, "%s.2.bt2" % org, "%s.3.bt2" % org, \
		"%s.4.bt2" % org, "%s.rev.1.bt2" % org, "%s.rev.2.bt2" % org]:	
			if not check_program_exists(config["bt2idx"], ff):
				sys.exit(1)
	
	#test software works
	if args["software"]:
		print "======================Testing Rscript...======================"
		os.system("%s/Rscript --version" % config["Rscriptbin"])
		print "======================Testing python...======================"
		os.system("%s/python --version" % config["pythonbin"])
		print "======================Testing trimmomatic...======================"
		os.system("java -jar %s/%s -version" % (config["trimmomaticbin"], config["trimmomaticjarfile"]))
		print "======================Testing bowtie2...======================"
		os.system("%s/bowtie2 --version" % config["bowtie2bin"])
		print "======================Testing samtools...======================"
		os.system("%s/samtools --version" % config["samtoolsbin"])
		print "======================Testing picard...======================"
		os.system("java -jar %s/%s -h" % (config["picardbin"], config["picardjarfile"]))
		print "======================Testing macs2...======================"
		os.system("pythonlib=`echo $PYTHONPATH | tr : \\n | grep -v %s | paste -s -d:` && unset $PYTHONPATH && export PYTHONPATH=$pythonlib:%s && %s/macs2 --version" % (config["macs2pythonlib"], config["macs2pythonlib"], config["macs2bin"]))
		print "======================Testing kseq...======================"
		os.system("%s/kseq_test --help" % config["kseqbin"])
		print "======================Testing meme...======================"
		os.system("%s/meme -version" % config["memebin"])
		print "======================Testing bedops...======================"
		os.system("%s/bedops --version" % config["bedopsbin"])
		print "======================Testing bedtools...======================"
		os.system("%s/bedtools --version" % config["bedtoolsbin"])
		print "======================Testing make_cut_matrix...======================"
		os.system("%s/make_cut_matrix --version" % config["makecutmatrixbin"])

