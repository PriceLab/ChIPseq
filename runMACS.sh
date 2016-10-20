# perform peak-calling with MACS
# Smad3


datadir=/proj/price1/jpearl/ChIPseq/
homedir=/proj/price1/sament/chipseq

# create files with Smad3 and control reads
# used as the input to MACS

echo merging Smad3 samples
cat $datadir/*/*Smad3*bed.sorted > $homedir/Smad3.merged.pairs.bed
cd $homedir
sort -k1,1 -k2,2n Smad3.merged.pairs.bed > Smad3.merged.pairs.sorted.bed

echo merging input controls
cat $datadir/*/*IC*bed.sorted > $homedir/IC.merged.pairs.bed
cd $homedir
sort -k1,1 -k2,2n $homedir/IC.merged.pairs.bed > IC.merged.pairs.sorted.bed

echo merging PolII samples
cat $datadir/*/*RNP2*bed.sorted > $homedir/RNP2.merged.pairs.bed
cd $homedir
sort -k1,1 -k2,2n RNP2.merged.pairs.bed > RNP2.merged.pairs.sorted.bed


# call peaks on pooled samples
# use lenient threshold as input to downstream
# filtering steps (e.g., IDR)

macs2 callpeak -t Smad3.merged.pairs.sorted.bed \
-c IC.merged.pairs.sorted.bed \
-f BEDPE -g mm -p .01 -n Smad3.combined.MACS.p01 \
--to-large --call-summits 

# call peaks on individual samples
# use lenient thresholds as input to downstream 
# filtering steps (e.g., IDR)

cp $datadir/*/*Smad3*.bed.sorted $homedir
for file in *Smad3*.bed.sorted
do
echo $file
macs2 callpeak -t $file \
-c IC.merged.pairs.sorted.bed \
-f BEDPE -g mm -p 1e-3 -n $file.MACS.p1e3 \
--to-large
done

# filter peaks with FDR > 0.01

for file in *p1e3*.narrowPeak
do
echo $file
awk 'BEGIN{OFS="\t"}{if ($9 > 2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }' \
$file > $file.fdr01.narrowPeak 
done

# merge the peaks from individual samples

cat *sorted.MACS.p1e3_peaks.narrowPeak > MACS.combined.narrowPeak
sort -k1,1 -k2,2n MACS.combined.narrowPeak > MACS.sorted.bed
bedtools merge -i MACS.sorted.bed > MACS.merged.bed

cat *sorted.MACS.p1e3_peaks.narrowPeak.fdr01.narrowPeak \
> MACS.combined.fdr01.narrowPeak
sort -k1,1 -k2,2n MACS.combined.fdr01.narrowPeak > MACS.fdr01.sorted.bed

# consider three different merging strategies
# merge peaks with no gaps, 100bp gaps, 200bp gaps

bedtools merge -i MACS.fdr01.sorted.bed > MACS.fdr01.merged.bed
bedtools merge -i MACS.fdr01.sorted.bed -d 100 > MACS.fdr01.merged.100bp.bed
bedtools merge -i MACS.fdr01.sorted.bed -d 200 > MACS.fdr01.merged.200bp.bed


# we can also consider a more conservative peak-calling strategy

echo repeating peak-calling with FDR<0.05 and --to-small

for file in *Smad3*.bed.sorted
do
echo $file
macs2 callpeak -t $file \
-c IC.merged.pairs.downsampled.sorted.bed \
-f BEDPE -g mm -p 0.05 -n $file.MACS.p05.downsampledIC. \
--call-summits
done

# merge the peaks from individual samples

cat *sorted.MACS.p05.downsampledIC._peaks.narrowPeak > MACS.combined.downsampledIC.tosmall.narrowPeak
sort -k1,1 -k2,2n MACS.combined.downsampledIC.tosmall.narrowPeak > MACS.downsampledIC.tosmall.sorted.bed
bedtools merge -i MACS.downsampledIC.tosmall.sorted.bed -d 200 > MACS.downsampledIC.tosmall.merged.d200.bed

# in the end, this strategy seemed to perform less well 
# than the top strategy at predicting peaks
# the top strategy is consistent with the IDR framework
# so I'll go with that.


