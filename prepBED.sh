FILE=$1
awk '{print $1 "\t" $3-1 "\t" $4 "\t" $6 "\t" $5 "\t" $2}' $FILE.pairs > $FILE.bed
sort -k1,1 -k2,2n $FILE.bed > $FILE.bed.sorted
