#This is a simple bash script to change the data structure of the MAnorm peaks output so that it works in the other R scripts
grep "unique_peak1" $1  > MAnorm_uniquePeaks1.txt
grep "unique_peak2" $1  > MAnorm_uniquePeaks2.txt

cut MAnorm_uniquePeaks1.txt -f 1-5  > MAnorm_uniquePeaks1_cut.txt
cut MAnorm_uniquePeaks2.txt -f 1-5 > MAnorm_uniquePeaks2_cut.txt

rm MAnorm_uniquePeaks1.txt
rm MAnorm_uniquePeaks2.txt
