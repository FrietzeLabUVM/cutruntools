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

#usage: ./integrated.motif.find.sh CR_NFYA_BCKO_Day4_r1_S9_aligned_reads_peaks.narrowPeak
#input is a narrowPeak file from MACS2

memebin=/home/qz64/meme/bin
bedopsbin=/home/qz64/bin
bedtoolsbin=/n/app/bedtools/2.26.0/bin
genome_sequence=/home/qz64/chrom.hg19/hg19.fa

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

>&2 echo "Get filtered peaks..."
$bedopsbin/bedops --range 150 -u $workdir/$dname/$summit > padded/$summit
$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed padded/$summit -fo padded.fa/$summitfa
./filter.py padded.fa/$summitfa 300 > repeat.region/$summit
$bedopsbin/bedops -n 1 $workdir/$dname/$summit repeat.region/$summit | $bedopsbin/sort-bed - > filtered/$summit
$bedopsbin/bedops -e 1 $workdir/$dname/$peak filtered/$summit > filtered/$peak

#motif discovery starts here
motif_dir=random.5000
msummit=$motif_dir/summits
mpadded=$motif_dir/padded
mpaddedfa=$motif_dir/padded.fa

for d in $motif_dir $msummit $mpadded $mpaddedfa; do
if [ ! -d $d ]; then
mkdir $d
fi
done

>&2 echo "Get randomized 5000 peaks..."
cat filtered/$peak | sort -t"	" -g -k8 -r | head -n 15000 | shuf | head -n 5000 | $bedopsbin/sort-bed - > $motif_dir/$peak
$bedopsbin/bedops -e 1 filtered/$summit $motif_dir/$peak > $msummit/$summit
$bedopsbin/bedops --range 150 -u $msummit/$summit > $mpadded/$summit

$bedtoolsbin/bedtools getfasta -fi $genome_sequence -bed $mpadded/$summit -fo $mpaddedfa/$summitfa

>&2 echo "Start MEME analysis for de novo motif finding..."
meme_outdir=$motif_dir/MEME_"$fname"_shuf
$memebin/meme-chip -oc $meme_outdir -dreme-m 20 -meme-nmotifs 20 $mpaddedfa/$summitfa

>&2 echo "Finished"
