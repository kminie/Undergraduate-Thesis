#This is a simple bash script to extract the correct columns from the peak and read files so that they can be inputted into MAnorm
if [ $# -ne 4 ]
then
  echo "Usage: `basename $0` peak1.bed peak2.bed read1.bed read2.bed"
  exit 
fi

cut $1 -f 1-3  > peak1.bed
cut $2 -f 1-3 > peak2.bed
cut $3 -f 1-3,6 > read1.bed
cut $4 -f 1-3,6 > read2.bed
