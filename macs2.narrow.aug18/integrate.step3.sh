#!/bin/bash
#SBATCH -n 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=32000                        # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bernardzhu@gmail.com   # Email to which notifications will be sent

#==============================================
# remember to change chrom.hg19 in line genome_sequence to mouse mm10 sequence. (simply replace it with /home/qz64/chrom.mm10/mm10.fa)
# (here I used the masked sequence file, masked means RepeatMasked)
#==============================================

#usage: ./integrated.step3.sh CR_NFYA_BCKO_Day4_r1_S9_aligned_reads_peaks.narrowPeak
#input is a narrowPeak file from MACS2

peak_file=$1 #a narrowPeak file
mbase=`basename $peak_file _peaks.narrowPeak`
mdiscovery=random.5000/MEME_"$mbase"_shuf

./read.meme.py $mdiscovery

memebin=/home/qz64/meme/bin
bedopsbin=/home/qz64/bin
bedtoolsbin=/n/app/bedtools/2.26.0/bin
genome_sequence=/home/qz64/chrom.hg19/hg19.fa

p=0.0005
motif_dir=$mdiscovery/motifs #a directory containing a list of *.meme files
base=`basename $peak_file .narrowPeak`
workdir=`pwd`
dir=`dirname $peak_file`
fa_dir=filtered.fa

if [ ! -d $fa_dir ]; then
mkdir $fa_dir
fi

$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $workdir/$dir/"$base".narrowPeak -fo $fa_dir/"$base".fa
./fix_sequence.py $fa_dir/"$base".fa

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

