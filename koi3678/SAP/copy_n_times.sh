#!/bin/bash
#

echo 'max'
read max

# Count number of files which need detrending
for i in `seq 2 $max`; do
  newpolyam=`echo polyam"$i".nb`
  newgp=`echo gp"$i".nb`
	newlocal=`echo local"$i".nb`
  cp polyam1.nb $newpolyam
	cp gp1.nb $newgp
	cp local1.nb $newlocal
done