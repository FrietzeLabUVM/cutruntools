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
javabin=%s

bt2idx=%s
kseqbin=%s

infile=$1
#expand the path of infile
relinfile=`realpath -s $infile`
dirname=`dirname $relinfile`
base=`basename $infile _R1_001.fastq.gz`
>&2 echo "Input file is $relinfile"
>&2 date

#cd to current directory
cd $dirname
workdir=`pwd`

len=`cat length`
trimdir=$workdir/trimmed
trimdir2=$workdir/trimmed3
logdir=$workdir/logs
aligndir=$workdir/aligned.aug10

for d in $trimdir $trimdir2 $logdir $aligndir; do
if [ ! -d $d ]; then
mkdir $d
fi
done

""" % (config["trimmomaticbin"], config["trimmomaticjarfile"], config["adapterpath"], 
	config["bowtie2bin"], config["samtoolsbin"], config["javabin"],
	config["bt2idx"], config["kseqbin"])

	bowtie2_org = config["input/output"]["organism_build"]
	if config["input/output"]["organism_build"]=="hg38":
		bowtie2_org = "GRCh38"

	scripts="""
#trimming paired-end
#good version
>&2 echo "Trimming file $base ..."
>&2 date
$javabin/java -jar $trimmomaticbin/$trimmomaticjarfile PE -threads 1 -phred33 $dirname/"$base"_R1_001.fastq.gz $dirname/"$base"_R2_001.fastq.gz $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

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

	p_pythonbase = config["pythonbin"].rstrip("/").rstrip("/bin")

	path = """picardbin=%s
picardjarfile=%s
samtoolsbin=%s
macs2bin=%s
javabin=%s
extratoolsbin=%s
extrasettings=%s
chromsizedir=`dirname %s`
macs2pythonlib=%s

pythonlib=`echo $PYTHONPATH | tr : "\\n" | grep -v $macs2pythonlib | paste -s -d:`
unset PYTHONPATH
export PYTHONPATH=$macs2pythonlib:$pythonlib

pythonldlibrary=%s
ldlibrary=`echo $LD_LIBRARY_PATH | tr : "\\n" | grep -v $pythonldlibrary | paste -s -d:`
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$pythonldlibrary:$ldlibrary

""" % (config["picardbin"], config["picardjarfile"], config["samtoolsbin"], 
	config["macs2bin"], config["javabin"],
	config["extratoolsbin"], config["extrasettings"], 
	config["genome_sequence"], config["macs2pythonlib"], p_pythonbase + "/lib")

	scripts =""">&2 echo "Input parameters are: $1"
>&2 date

#expand the path of $1
relinfile=`realpath -s $1`
dirname=`dirname $relinfile`
base=`basename $1 .bam`

#cd to current directory (aligned.aug10)
cd $dirname

workdir=`pwd`
logdir=$workdir/logs

for d in $logdir sorted dup.marked dedup; do
if [ ! -d $d ]; then
mkdir $d
fi
done

>&2 echo "Sorting bam... ""$base".bam
>&2 date
$javabin/java -jar $picardbin/$picardjarfile SortSam \
INPUT=$dirname/"$base".bam OUTPUT=sorted/"$base".bam SORT_ORDER=coordinate

>&2 echo "Marking duplicates... ""$base".bam
>&2 date
$javabin/java -jar $picardbin/$picardjarfile MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dup.marked/"$base".bam \
METRICS_FILE=metrics."$base".txt

>&2 echo "Removing duplicates... ""$base".bam
>&2 date
$javabin/java -jar $picardbin/$picardjarfile MarkDuplicates \
INPUT=sorted/"$base".bam OUTPUT=dedup/"$base".bam \
METRICS_FILE=metrics."$base".txt \
REMOVE_DUPLICATES=true

for d in sorted.120bp dup.marked.120bp dedup.120bp; do
if [ ! -d $d ]; then
mkdir $d
fi
done

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
perlbin=%s
genome_sequence=%s
extrasettings=%s
blacklist=$extrasettings/%s.blacklist.bed

i=$1 #filename must end with .narrowPeak
>&2 echo "Input file is $i"

#expand the path for $1
relinfile=`realpath -s $i`
dirname=`dirname $relinfile`

#cd to current directory (macs2.narrow.aug10)
cd $dirname

for d in padded padded.fa repeat.region filtered blk_filtered; do
if [ ! -d $d ]; then
mkdir $d
fi
done

workdir=`pwd`
fname=`basename $i _peaks.narrowPeak`
peak=$fname"_peaks.narrowPeak"
summit=$fname"_summits.bed"
summitfa=$fname"_summits_padded.fa"
""" % (config["memebin"], config["bedopsbin"], config["bedtoolsbin"], config["pythonbin"], config["perlbin"], config["genome_sequence"], config["extrasettings"], config["input/output"]["organism_build"])

	scripts2 = """>&2 echo "Get filtered peaks..."
$bedopsbin/bedops --range %d -u $workdir/$summit > padded/$summit
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed padded/$summit -fo padded.fa/$summitfa
$pythonbin/python filter.py padded.fa/$summitfa %d > repeat.region/$summit
$bedopsbin/bedops -n 1 $workdir/$summit repeat.region/$summit | $bedopsbin/sort-bed - > filtered/$summit
$bedopsbin/bedops -e 1 $workdir/$peak filtered/$summit > filtered/$peak
cat filtered/$peak | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > blk_filtered/$peak
cat filtered/$summit | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > blk_filtered/$summit

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
cat blk_filtered/$peak | sort -t"	" -g -k8 -r | head -n %d | shuf | head -n %d | $bedopsbin/sort-bed - > $motif_dir/$peak
$bedopsbin/bedops -e 1 blk_filtered/$summit $motif_dir/$peak > $msummit/$summit
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

#expand the path for $peak_file
relinfile=`realpath -s $peak_file`
dirname=`dirname $relinfile`

#cd to current directory (macs2.narrow.aug10)
cd $dirname


$pythonbin/python read.meme.py $mdiscovery
""" % (config["pythonbin"], config["motif_finding"]["num_peaks"])

	p_pythonbase = config["pythonbin"].rstrip("/").rstrip("/bin")

	scripts2 = """
memebin=%s
bedopsbin=%s
bedtoolsbin=%s
genome_sequence=%s
samtoolsbin=%s
makecutmatrixbin=%s
Rscriptbin=%s
extrasettings=%s

pythonldlibrary=%s
ldlibrary=`echo $LD_LIBRARY_PATH | tr : "\\n" | grep -v $pythonldlibrary | paste -s -d:`
unset LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$pythonldlibrary:$ldlibrary

p=%.5f
motif_dir=$mdiscovery/motifs #a directory containing a list of *.meme files
base=`basename $peak_file .narrowPeak`
workdir=`pwd`
dir=blk_filtered
fa_dir=blk_filtered.fa

if [ ! -d $fa_dir ]; then
mkdir $fa_dir
fi

blacklist=$extrasettings/%s.blacklist.bed
cat $workdir/$dir/"$base".narrowPeak | grep -v -e "chrM" | $bedopsbin/sort-bed - | $bedopsbin/bedops -n 1 - $blacklist > $workdir/$dir/"$base".filtered.narrowPeak
""" % (config["memebin"], config["bedopsbin"], config["bedtoolsbin"], 
	config["genome_sequence"], config["samtoolsbin"], 
	config["makecutmatrixbin"], config["Rscriptbin"],
	config["extrasettings"], p_pythonbase + "/lib", 
	config["motif_finding"]["motif_scanning_pval"], 
	config["input/output"]["organism_build"])

	scripts3 = """
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $workdir/$dir/"$base".filtered.narrowPeak -fo $fa_dir/"$base".fa
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

cur_path=`echo $PATH | tr : "\\n" | grep -v $bedopsbin | paste -s -d:`
unset PATH
export PATH=$cur_path:$bedopsbin

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

peakfile=blk_filtered/"$base"_peaks.narrowPeak
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

def generate_integrated_all_steps_sh(output=None):
	script = """#!/bin/bash

sample=$1
#expand the path for $sample
relsample=`realpath -s $sample`
dirname=`dirname $relsample`
base=`basename $relsample _R1_001.fastq.gz`

cd $dirname

echo "Submitting job 1..."
jid1=$(sbatch ./integrated.sh $relsample)
jid1=${jid1##* }
echo "Submitting job 2..."
jid2=$(sbatch -d afterany:$jid1 aligned.aug10/integrated.step2.sh aligned.aug10/"$base"_aligned_reads.bam)
jid2=${jid2##* }
echo "Submitting job 3..."
jid3=$(sbatch -d afterany:$jid2 macs2.narrow.aug18/integrate.motif.find.sh macs2.narrow.aug18/"$base"_aligned_reads_peaks.narrowPeak)
jid3=${jid3##* }
echo "Submitting job 4..."
jid4=$(sbatch -d afterany:$jid3 macs2.narrow.aug18/integrate.footprinting.sh macs2.narrow.aug18/"$base"_aligned_reads_peaks.narrowPeak)
jid4=${jid4##* }

echo "Submitting job 5..."
jid5=$(sbatch -d afterany:$jid4 macs2.narrow.aug18.dedup/integrate.motif.find.sh macs2.narrow.aug18.dedup/"$base"_aligned_reads_peaks.narrowPeak)
jid5=${jid5##* }
echo "Submitting job 6..."
jid6=$(sbatch -d afterany:$jid5 macs2.narrow.aug18.dedup/integrate.footprinting.sh macs2.narrow.aug18.dedup/"$base"_aligned_reads_peaks.narrowPeak)
jid6=${jid6##* }
"""
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw
	outp.write(script + "\n")
	if output is not None:
		outp.close()


def generate_single_locus_script_sh(config, dedup=False, output=None):
	outp = sys.stdout
	if output is not None:
		fw = open(output, "w")
		outp = fw

	if dedup:
		samtools_flags = "-f 3 -F 4 -F 8 -F 1024"
	else:
		samtools_flags = "-f 3 -F 4 -F 8"

	script = """#!/bin/bash

region=$1
bamfile=$2
outdir=$3
chromsizedir=`dirname %s`
chromsizefile=$chromsizedir/%s.chrom.sizes
pythonbin=%s
samtoolsbin=%s
bedtoolsbin=%s
bedopsbin=%s
extratoolsbin=%s
samtoolsflags="%s"

regionname=`echo $region|sed "s/:/-/g"`
basename=`basename $bamfile .bam`
newbamfile="$basename"-"$regionname".bam
newbase=`basename $newbamfile .bam`

if [ ! -d $outdir ]; then
mkdir $outdir
fi
$samtoolsbin/samtools view -bh $samtoolsflags $bamfile "$region" > $outdir/$newbamfile
$samtoolsbin/samtools index $outdir/$newbamfile
$samtoolsbin/samtools view -b $outdir/$newbamfile|$samtoolsbin/samtools sort -O bam -n - -T tmp.test|$bedtoolsbin/bedtools bamtobed -i stdin -bedpe > $outdir/"$newbase".frag.ends.txt

$pythonbin/python check_coordinate.py $chromsizefile $outdir/"$newbase".frag.ends.txt > $outdir/"$newbase".frag.ends.checked.txt

$pythonbin/python quantify_separate.py $outdir/"$newbase".frag.ends.checked.txt $outdir/"$newbase".frag.ends.R1.bed $outdir/"$newbase".frag.ends.R2.bed
$bedopsbin/sort-bed $outdir/"$newbase".frag.ends.R1.bed > $outdir/"$newbase".frag.ends.R1.sorted.bed
$bedopsbin/sort-bed $outdir/"$newbase".frag.ends.R2.bed > $outdir/"$newbase".frag.ends.R2.sorted.bed
$bedtoolsbin/groupBy -i $outdir/"$newbase".frag.ends.R1.sorted.bed -g 1,2,3 -c 2 -o count > $outdir/"$newbase".frag.ends.R1.bdg
$bedtoolsbin/groupBy -i $outdir/"$newbase".frag.ends.R2.sorted.bed -g 1,2,3 -c 2 -o count > $outdir/"$newbase".frag.ends.R2.bdg
$extratoolsbin/bedGraphToBigWig $outdir/"$newbase".frag.ends.R1.bdg $chromsizefile $outdir/"$newbase".frag.ends.R1.bw
$extratoolsbin/bedGraphToBigWig $outdir/"$newbase".frag.ends.R2.bdg $chromsizefile $outdir/"$newbase".frag.ends.R2.bw

$pythonbin/python quantify.py $outdir/"$newbase".frag.ends.checked.txt $outdir/"$newbase".frag.ends.bed
$bedopsbin/sort-bed $outdir/"$newbase".frag.ends.bed > $outdir/"$newbase".frag.ends.sorted.bed
$bedtoolsbin/groupBy -i $outdir/"$newbase".frag.ends.sorted.bed -g 1,2,3 -c 2 -o count > $outdir/"$newbase".frag.ends.bdg
$extratoolsbin/bedGraphToBigWig $outdir/"$newbase".frag.ends.bdg $chromsizefile $outdir/"$newbase".frag.ends.bw

""" % (config["genome_sequence"], config["input/output"]["organism_build"], config["pythonbin"], config["samtoolsbin"], config["bedtoolsbin"], \
config["bedopsbin"], config["extratoolsbin"], samtools_flags)

	outp.write(script + "\n")
	if output is not None:
		outp.close()



def write_length_file(n, length):
	fw = open(n, "w")
	fw.write(str(length) + "\n")
	fw.close()

if __name__=="__main__":
	curpath = os.path.dirname(os.path.abspath(__file__))
	#print curpath
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
	generate_integrated_all_steps_sh(output=outdir+"/integrated.all.steps.sh")
	generate_single_locus_script_sh(config, dedup=False, output=outdir+"/macs2.narrow.aug18/get_cuts_single_locus.sh")
	generate_single_locus_script_sh(config, dedup=True, output=outdir+"/macs2.narrow.aug18.dedup/get_cuts_single_locus.sh")

	make_executable(outdir+"/integrated.sh")
	make_executable(outdir+"/integrated.all.steps.sh")
	make_executable(outdir+"/aligned.aug10/integrated.step2.sh")
	make_executable(outdir+"/macs2.narrow.aug18/integrate.motif.find.sh")
	make_executable(outdir+"/macs2.narrow.aug18.dedup/integrate.motif.find.sh")
	make_executable(outdir+"/macs2.narrow.aug18/integrate.footprinting.sh")
	make_executable(outdir+"/macs2.narrow.aug18.dedup/integrate.footprinting.sh")
	make_executable(outdir+"/macs2.narrow.aug18/get_cuts_single_locus.sh")
	make_executable(outdir+"/macs2.narrow.aug18.dedup/get_cuts_single_locus.sh")

	shutil.copyfile("%s/macs2.narrow.aug18/filter.py" % curpath, outdir+"/macs2.narrow.aug18/filter.py")
	shutil.copyfile("%s/macs2.narrow.aug18/fix_sequence.py" % curpath, outdir+"/macs2.narrow.aug18/fix_sequence.py")
	shutil.copyfile("%s/macs2.narrow.aug18/read.meme.py" % curpath, outdir+"/macs2.narrow.aug18/read.meme.py")
	shutil.copyfile("%s/macs2.narrow.aug18/run_centipede_parker.R" % curpath, outdir+"/macs2.narrow.aug18/run_centipede_parker.R")
	shutil.copyfile("%s/quantify.py" % curpath, outdir+"/macs2.narrow.aug18/quantify.py")
	shutil.copyfile("%s/quantify_separate.py" % curpath, outdir+"/macs2.narrow.aug18/quantify_separate.py")
	shutil.copyfile("%s/check_coordinate.py" % curpath, outdir+"/macs2.narrow.aug18/check_coordinate.py")

	shutil.copyfile("%s/macs2.narrow.aug18/filter.py" % curpath, outdir+"/macs2.narrow.aug18.dedup/filter.py")
	shutil.copyfile("%s/macs2.narrow.aug18/fix_sequence.py" % curpath, outdir+"/macs2.narrow.aug18.dedup/fix_sequence.py")
	shutil.copyfile("%s/macs2.narrow.aug18/read.meme.py" % curpath, outdir+"/macs2.narrow.aug18.dedup/read.meme.py")
	shutil.copyfile("%s/macs2.narrow.aug18/run_centipede_parker.R" % curpath, outdir+"/macs2.narrow.aug18.dedup/run_centipede_parker.R")
	shutil.copyfile("%s/quantify.py" % curpath, outdir+"/macs2.narrow.aug18.dedup/quantify.py")
	shutil.copyfile("%s/quantify_separate.py" % curpath, outdir+"/macs2.narrow.aug18.dedup/quantify_separate.py")
	shutil.copyfile("%s/check_coordinate.py" % curpath, outdir+"/macs2.narrow.aug18.dedup/check_coordinate.py")

	shutil.copyfile(sys.argv[1], outdir+"/config.json")

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
