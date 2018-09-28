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

#===============================================================================================
#Usage: ./integrated.step4.sh ../aligned.aug10/dup.marked.120bp/CR_NFYA_RNP115_Day4_r1_S12_aligned_reads.bam
#(input file please use the BAM file from ../aligned.aug10/dup.marked.120bp/*.bam, ie the dup.marked.120bp version
#===============================================================================================

module load samtools

bamfile=$1 #use the dup.marked.120bp version
workdir=`pwd`
dir=`dirname $bamfile`
base=`basename $bamfile .bam`

dest=centipede.bam
outbam=$dest/"$base".bam
if [ ! -d $dest ]; then
mkdir $dest
fi
#note that 1024 means read is PCR or optical duplicate
#samtools view -b -h -f 3 -F 4 -F 8 -F 1024 -o $outbam $bamfile #previous version
samtools view -b -h -f 3 -F 4 -F 8 -o $outbam $bamfile
samtools index $outbam

peakfile=filtered/"$base"_peaks.narrowPeak
fimo_dir=fimo.result/"$base"_peaks

for i in `ls -1 $fimo_dir`; do #shows a list of motifs
echo "Doing $i..."
fimo_d=$fimo_dir/$i
tmp=`echo $i|cut -d "." -f3|wc -c`
mlen=$(( tmp - 1 ))
~/.local/bin/make_cut_matrix -v -b '(25-150 1)' -d -o 0 -r 100 -p 1 -f 3 -F 4 -F 8 -q 0 $outbam $fimo_d/fimo.bed > $fimo_d/fimo.cuts.freq.txt
./run_centipede_parker.R $fimo_d/fimo.cuts.freq.txt $fimo_d/fimo.bed $fimo_d/fimo.png $mlen
done

