#!/usr/bin/python
import shutil
import sys
import os
import re
import json

def make_executable(filename):
	st = os.stat(filename)
	os.chmod(filename, st.st_mode | 0o111)
	
def generate_integrated_sh(config, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
		
	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%d                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_alignment"]["time_limit"], config["cluster"]["step_alignment"]["queue"], 
	config["cluster"]["step_alignment"]["memory"], config["cluster"]["email"])

	path="""trimmomaticbin=%s
trimmomaticjarfile=%s
adapterpath=%s
bowtie2bin=%s
samtoolsbin=%s
len=`cat length`

workdir=`pwd`
bt2idx=%s
trimdir=$workdir/trimmed
trimdir2=$workdir/trimmed3
logdir=$workdir/logs
aligndir=$workdir/aligned.aug10
kseqbin=%s

mkdir $trimdir
mkdir $trimdir2
mkdir $logdir
mkdir $aligndir

infile=$1
dirname=`dirname $infile`
base=`basename $infile _R1_001.fastq.gz`
""" % (config["trimmomaticbin"], config["trimmomaticjarfile"], config["adapterpath"], 
	config["bowtie2bin"], config["samtoolsbin"], 
	config["bt2idx"], config["kseqbin"])

	bowtie2_org = config["input/output"]["organism_build"]
	if config["input/output"]["organism_build"]=="hg38":
		bowtie2_org = "GRCh38"

	scripts=""">&2 echo "Input file is $infile"
>&2 date

#trimming paired-end
#good version
>&2 echo "Trimming file $base ..."
>&2 date
java -jar $trimmomaticbin/$trimmomaticjarfile PE -threads 1 -phred33 $dirname/"$base"_R1_001.fastq.gz $dirname/"$base"_R2_001.fastq.gz $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

>&2 echo "Second stage trimming $base ..."
>&2 date
$kseqbin/kseq_test $trimdir/"$base"_1.paired.fastq.gz $len $trimdir2/"$base"_1.paired.fastq.gz
$kseqbin/kseq_test $trimdir/"$base"_2.paired.fastq.gz $len $trimdir2/"$base"_2.paired.fastq.gz

>&2 echo "Aligning file $base ..."
>&2 date
($bowtie2bin/bowtie2 -p 2 --dovetail --phred33 -x $bt2idx/%s -1 $trimdir2/"$base"_1.paired.fastq.gz -2 $trimdir2/"$base"_2.paired.fastq.gz) 2> $logdir/"$base".bowtie2 | $samtoolsbin/samtools view -bS - > $aligndir/"$base"_aligned_reads.bam

>&2 echo "Finished"
>&2 date
""" % bowtie2_org
	outp.write(header + "\n")
	outp.write(path + "\n")
	outp.write(scripts + "\n")

	if output is not None:
		outp.close()

def generate_integrated_step2_sh(config, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw

	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%d                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_process_bam"]["time_limit"], config["cluster"]["step_process_bam"]["queue"], 
	config["cluster"]["step_process_bam"]["memory"], config["cluster"]["email"])

	path = """picardbin=%s
picardjarfile=%s
samtoolsbin=%s
macs2bin=%s
extratoolsbin=%s
extrasettings=%s
chromsizedir=`dirname %s`
macs2pythonlib=%s
workdir=`pwd`
logdir=$workdir/logs

mkdir $logdir
pythonlib=`echo $PYTHONPATH | tr : "\\n" | grep -v $macs2pythonlib | paste -s -d:`
unset $PYTHONPATH
export PYTHONPATH=$pythonlib:$macs2pythonlib

""" % (config["picardbin"], config["picardjarfile"], config["samtoolsbin"], 
	config["macs2bin"], config["extratoolsbin"], config["extrasettings"], 
	config["genome_sequence"], config["macs2pythonlib"])

	scripts =""">&2 echo "Input parameters are: $1"
>&2 date
base=`basename $1 .bam`

mkdir sorted dup.marked dedup
>&2 echo "Sorting bam... ""$base".bam
>&2 date
java -jar $picardbin/$picardjarfile SortSam \
INPUT="$base".bam OUTPUT=sorted/"$base".bam SORT_ORDER=coordinate

>&2 echo "Marking duplicates... ""$base".bam
>&2 date
java -jar $picardbin/$picardjarfile MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dup.marked/"$base".bam \
METRICS_FILE=metrics."$base".txt

>&2 echo "Removing duplicates... ""$base".bam
>&2 date
java -jar $picardbin/$picardjarfile MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dedup/"$base".bam \
METRICS_FILE=metrics."$base".txt \
REMOVE_DUPLICATES=true

mkdir sorted.120bp dup.marked.120bp dedup.120bp
>&2 echo "Filtering to <120bp... ""$base".bam
>&2 date
$samtoolsbin/samtools view -h sorted/"$base".bam |awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > sorted.120bp/"$base".bam
$samtoolsbin/samtools view -h dup.marked/"$base".bam |awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dup.marked.120bp/"$base".bam
$samtoolsbin/samtools view -h dedup/"$base".bam |awk -f $extrasettings/filter_below.awk |$samtoolsbin/samtools view -Sb - > dedup.120bp/"$base".bam

>&2 echo "Creating bam index files... ""$base".bam
>&2 date
$samtoolsbin/samtools index sorted/"$base".bam
$samtoolsbin/samtools index dup.marked/"$base".bam
$samtoolsbin/samtools index dedup/"$base".bam
$samtoolsbin/samtools index sorted.120bp/"$base".bam
$samtoolsbin/samtools index dup.marked.120bp/"$base".bam
$samtoolsbin/samtools index dedup.120bp/"$base".bam

>&2 echo "Peak calling using MACS2... ""$base".bam
>&2 echo "Logs are stored in $logdir"
>&2 date
bam_file=dup.marked.120bp/"$base".bam
dir=`dirname $bam_file`
base_file=`basename $bam_file .bam`
"""

	macs2_org = "hs"
	if config["input/output"]["organism_build"]=="hg19" or \
	config["input/output"]["organism_build"]=="hg38":
		macs2_org = "hs"
	elif config["input/output"]["organism_build"]=="mm10" or \
	config["input/output"]["organism_build"]=="mm9":
		macs2_org = "mm"

	macs_script = """
outdir=$workdir/../macs2.narrow.aug18 #for macs2
outdir2=$workdir/../macs2.narrow.aug18.dedup #for macs2 dedup version

for d in $outdir $outdir2; do
if [ ! -d $d ]; then
mkdir $d
fi
done

$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdir -q 0.01 -B --SPMR --keep-dup all 2> $logdir/"$base_file".macs2
$macs2bin/macs2 callpeak -t $workdir/$dir/"$base_file".bam -g %s -f BAMPE -n $base_file --outdir $outdir2 -q 0.01 -B --SPMR 2> $logdir/"$base_file".dedup.macs2
""" % (macs2_org, macs2_org)

	scripts2 = """
>&2 echo "Converting bedgraph to bigwig... ""$base".bam
>&2 date
cd $outdir
sort -k1,1 -k2,2n $outdir/"$base_file"_treat_pileup.bdg > $outdir/"$base_file".sort.bdg
$extratoolsbin/bedGraphToBigWig $outdir/"$base_file".sort.bdg $chromsizedir/%s.chrom.sizes $outdir/"$base_file".sorted.bw
rm -rf "$base_file".sort.bdg

cd $outdir2
sort -k1,1 -k2,2n $outdir2/"$base_file"_treat_pileup.bdg > $outdir2/"$base_file".sort.bdg
$extratoolsbin/bedGraphToBigWig $outdir2/"$base_file".sort.bdg $chromsizedir/%s.chrom.sizes $outdir2/"$base_file".sorted.bw
rm -rf "$base_file".sort.bdg

>&2 echo "Finished"
>&2 date
""" % (config["input/output"]["organism_build"], config["input/output"]["organism_build"])

	outp.write(header + "\n")
	outp.write(path + "\n")
	outp.write(scripts + "\n")
	outp.write(macs_script + "\n")
	outp.write(scripts2 + "\n")

	if output is not None:
		outp.close()

def generate_integrated_motif_find_sh(config, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
	
	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%s                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_motif_find"]["time_limit"], config["cluster"]["step_motif_find"]["queue"], 
	config["cluster"]["step_motif_find"]["memory"], config["cluster"]["email"])
	
	scripts = """memebin=%s
bedopsbin=%s
bedtoolsbin=%s
pythonbin=%s
genome_sequence=%s

i=$1 #filename must end with .narrowPeak
>&2 echo "Input file is $i"

for d in padded padded.fa repeat.region filtered; do
if [ ! -d $d ]; then
mkdir $d
fi
done

workdir=`pwd`

dname=`dirname $i`
fname=`basename $i _peaks.narrowPeak`
peak=$fname"_peaks.narrowPeak"
summit=$fname"_summits.bed"
summitfa=$fname"_summits_padded.fa"
""" % (config["memebin"], config["bedopsbin"], config["bedtoolsbin"], config["pythonbin"], config["genome_sequence"])


	scripts2 = """>&2 echo "Get filtered peaks..."
$bedopsbin/bedops --range %d -u $workdir/$dname/$summit > padded/$summit
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed padded/$summit -fo padded.fa/$summitfa
$pythonbin/python filter.py padded.fa/$summitfa %d > repeat.region/$summit
$bedopsbin/bedops -n 1 $workdir/$dname/$summit repeat.region/$summit | $bedopsbin/sort-bed - > filtered/$summit
$bedopsbin/bedops -e 1 $workdir/$dname/$peak filtered/$summit > filtered/$peak
""" % (config["motif_finding"]["num_bp_from_summit"], 2*config["motif_finding"]["num_bp_from_summit"])

	scripts3 = """#motif discovery starts here
motif_dir=random.%d
msummit=$motif_dir/summits
mpadded=$motif_dir/padded
mpaddedfa=$motif_dir/padded.fa

for d in $motif_dir $msummit $mpadded $mpaddedfa; do
if [ ! -d $d ]; then
mkdir $d
fi
done

>&2 echo "Get randomized %d peaks..."
cat filtered/$peak | sort -t"	" -g -k8 -r | head -n %d | shuf | head -n %d | $bedopsbin/sort-bed - > $motif_dir/$peak
$bedopsbin/bedops -e 1 filtered/$summit $motif_dir/$peak > $msummit/$summit
$bedopsbin/bedops --range %d -u $msummit/$summit > $mpadded/$summit

$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $mpadded/$summit -fo $mpaddedfa/$summitfa

>&2 echo "Start MEME analysis for de novo motif finding..."
meme_outdir=$motif_dir/MEME_"$fname"_shuf
$memebin/meme-chip -oc $meme_outdir -dreme-m %d -meme-nmotifs %d $mpaddedfa/$summitfa

>&2 echo "Finished"
""" % (config["motif_finding"]["num_peaks"], config["motif_finding"]["num_peaks"], 
	config["motif_finding"]["total_peaks"], 
	config["motif_finding"]["num_peaks"], config["motif_finding"]["num_bp_from_summit"], 
	config["motif_finding"]["num_motifs"], config["motif_finding"]["num_motifs"])

	outp.write(header + "\n")
	outp.write(scripts + "\n")
	outp.write(scripts2 + "\n")
	outp.write(scripts3 + "\n")

	if output is not None:
		outp.close()

def generate_integrated_footprinting_sh(config, dedup=False, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw

	header = """#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t %s                         # Runtime in D-HH:MM format
#SBATCH -p %s                           # Partition to run in
#SBATCH --mem=%s                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=%s   # Email to which notifications will be sent
""" % (config["cluster"]["step_footprinting"]["time_limit"], config["cluster"]["step_footprinting"]["queue"], 
	config["cluster"]["step_footprinting"]["memory"], config["cluster"]["email"])

	scripts = """
pythonbin=%s
peak_file=$1 #a narrowPeak file
mbase=`basename $peak_file _peaks.narrowPeak`
mdiscovery=random.%d/MEME_"$mbase"_shuf

$pythonbin/python read.meme.py $mdiscovery
""" % (config["pythonbin"], config["motif_finding"]["num_peaks"])

	scripts2 = """
memebin=%s
bedopsbin=%s
bedtoolsbin=%s
genome_sequence=%s
samtoolsbin=%s
makecutmatrixbin=%s
Rscriptbin=%s

p=%.5f
motif_dir=$mdiscovery/motifs #a directory containing a list of *.meme files
base=`basename $peak_file .narrowPeak`
workdir=`pwd`
dir=`dirname $peak_file`
fa_dir=filtered.fa

if [ ! -d $fa_dir ]; then
mkdir $fa_dir
fi
""" % (config["memebin"], config["bedopsbin"], config["bedtoolsbin"], 
	config["genome_sequence"], config["samtoolsbin"], 
	config["makecutmatrixbin"], config["Rscriptbin"],
	config["motif_finding"]["motif_scanning_pval"])

	scripts3 = """
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $workdir/$dir/"$base".narrowPeak -fo $fa_dir/"$base".fa
$pythonbin/python fix_sequence.py $fa_dir/"$base".fa

outdir=fimo.result
for d in $outdir $outdir/$base; do
if [ ! -d $d ]; then
mkdir $d
fi
done

for m in `ls -1 $motif_dir`; do
motif=`basename $m .meme`
fimo_d=$outdir/$base/fimo2.$motif
if [ ! -d $fimo_d ]; then
mkdir $fimo_d
fi
$memebin/fimo --thresh $p --parse-genomic-coord -oc $fimo_d $motif_dir/"$motif".meme $fa_dir/"$base".fa
$bedopsbin/gff2bed < $fimo_d/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $fimo_d/fimo.bed
done
"""
	scripts4 = """
bamfile=../aligned.aug10/dup.marked.120bp/"$mbase".bam
workdir=`pwd`
dir=`dirname $bamfile`
base=`basename $bamfile .bam`

dest=centipede.bam
outbam=$dest/"$base".bam
if [ ! -d $dest ]; then
mkdir $dest
fi
"""
	scripts5 = """
$samtoolsbin/samtools view -b -h -f 3 -F 4 -F 8 -o $outbam $bamfile
"""
	if dedup:
		scripts5 = """
#note that 1024 means read is PCR or optical duplicate
$samtoolsbin/samtools view -b -h -f 3 -F 4 -F 8 -F 1024 -o $outbam $bamfile #previous version
"""

	scripts6 = """
$samtoolsbin/samtools index $outbam

peakfile=filtered/"$base"_peaks.narrowPeak
fimo_dir=fimo.result/"$base"_peaks

for i in `ls -1 $fimo_dir`; do #shows a list of motifs
echo "Doing $i..."
fimo_d=$fimo_dir/$i
tmp=`echo $i|cut -d "." -f3|wc -c`
mlen=$(( tmp - 1 ))
$makecutmatrixbin/make_cut_matrix -v -b '(25-150 1)' -d -o 0 -r 100 -p 1 -f 3 -F 4 -F 8 -q 0 $outbam $fimo_d/fimo.bed > $fimo_d/fimo.cuts.freq.txt
$Rscriptbin/Rscript run_centipede_parker.R $fimo_d/fimo.cuts.freq.txt $fimo_d/fimo.bed $fimo_d/fimo.png $mlen
done
"""
	outp.write(header + "\n")
	outp.write(scripts + "\n")
	outp.write(scripts2 + "\n")
	outp.write(scripts3 + "\n")
	outp.write(scripts4 + "\n")
	outp.write(scripts5 + "\n")
	outp.write(scripts6 + "\n")

	if output is not None:
		outp.close()

def write_length_file(n, length):
	fw = open(n, "w")
	fw.write(str(length) + "\n")
	fw.close()

if __name__=="__main__":
	f = open(sys.argv[1]) #config.json
	config = json.load(f)
	f.close()

	outdir = config["input/output"]["workdir"]
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	write_length_file(config["input/output"]["workdir"] + "/length", config["input/output"]["fastq_sequence_length"])

	if not os.path.isdir(outdir + "/aligned.aug10"):
		os.mkdir(outdir + "/aligned.aug10")
	if not os.path.isdir(outdir + "/macs2.narrow.aug18"):
		os.mkdir(outdir + "/macs2.narrow.aug18")
	if not os.path.isdir(outdir + "/macs2.narrow.aug18.dedup"):
		os.mkdir(outdir + "/macs2.narrow.aug18.dedup")

	generate_integrated_sh(config, output=outdir+"/integrated.sh")
	generate_integrated_step2_sh(config, output=outdir+"/aligned.aug10/integrated.step2.sh")
	generate_integrated_motif_find_sh(config, output=outdir+"/macs2.narrow.aug18/integrate.motif.find.sh")
	generate_integrated_motif_find_sh(config, output=outdir+"/macs2.narrow.aug18.dedup/integrate.motif.find.sh")
	generate_integrated_footprinting_sh(config, dedup=False, output=outdir+"/macs2.narrow.aug18/integrate.footprinting.sh")
	generate_integrated_footprinting_sh(config, dedup=True, output=outdir+"/macs2.narrow.aug18.dedup/integrate.footprinting.sh")

	make_executable(outdir+"/integrated.sh")
	make_executable(outdir+"/aligned.aug10/integrated.step2.sh")
	make_executable(outdir+"/macs2.narrow.aug18/integrate.motif.find.sh")
	make_executable(outdir+"/macs2.narrow.aug18.dedup/integrate.motif.find.sh")
	make_executable(outdir+"/macs2.narrow.aug18/integrate.footprinting.sh")
	make_executable(outdir+"/macs2.narrow.aug18.dedup/integrate.footprinting.sh")

	shutil.copyfile("macs2.narrow.aug18/filter.py", outdir+"/macs2.narrow.aug18/filter.py")
	shutil.copyfile("macs2.narrow.aug18/fix_sequence.py", outdir+"/macs2.narrow.aug18/fix_sequence.py")
	shutil.copyfile("macs2.narrow.aug18/read.meme.py", outdir+"/macs2.narrow.aug18/read.meme.py")
	shutil.copyfile("macs2.narrow.aug18/run_centipede_parker.R", outdir+"/macs2.narrow.aug18/run_centipede_parker.R")

	shutil.copyfile("macs2.narrow.aug18/filter.py", outdir+"/macs2.narrow.aug18.dedup/filter.py")
	shutil.copyfile("macs2.narrow.aug18/fix_sequence.py", outdir+"/macs2.narrow.aug18.dedup/fix_sequence.py")
	shutil.copyfile("macs2.narrow.aug18/read.meme.py", outdir+"/macs2.narrow.aug18.dedup/read.meme.py")
	shutil.copyfile("macs2.narrow.aug18/run_centipede_parker.R", outdir+"/macs2.narrow.aug18.dedup/run_centipede_parker.R")


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

	if config["input/output"]["fastq_directory"]!=config["input/output"]["workdir"]:
		for fx in flist:
			if os.path.islink(config["input/output"]["workdir"] + "/" + fx):
				os.remove(config["input/output"]["workdir"] + "/" + fx)
			os.symlink(config["input/output"]["fastq_directory"] + "/" + fx, \
			config["input/output"]["workdir"] + "/" + fx)
