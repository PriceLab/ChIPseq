#!/bin/bash

homedir=/proj/price1/sament/chipseq/Smad3
cd $homedir

bedfiles=/proj/price1/jpearl/ChIPseq/*/*Smad3*bed.sorted

wc -l $bedfiles > libsize.txt

meansize="$(awk 'NR<=4 { sum+=$1/4 } END {print int(sum) }' libsize.txt)"

chromsizes=/proj/price1/sament/resources/refGenomes/mm9/mm9.chrom.sizes

for i in {1..4}
do
  file="$(head -n $i libsize.txt | tail -n 1 | awk 'BEGIN{FS=" "} {print $2}' )"
  size="$(head -n $i libsize.txt | tail -n 1 | awk 'BEGIN{FS=" "} {print $1}' )"
  sizefactor="$(echo "scale = 10; $meansize / $size" | bc)"
  echo $file
  echo $sizefactor
  bedtools genomecov -bg -i $file -scale $sizefactor -g $chromsizes > $file.bedGraph
  bedGraphToBigWig $file.bedGraph $chromsizes $file.bw
  mv $file.bedGraph $file.bw $homedir
done
