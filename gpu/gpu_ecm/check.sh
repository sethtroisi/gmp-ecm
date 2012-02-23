#!/bin/sh

if [ -z $2 ] ; then
  n=0
else
  n=$2
fi

if [ -z $1 ] ; then
  B1=1024
else
  B1=$1
fi

if [ -z $3 ] ; then
  s=42424242
else
  s=$3
fi

N=32988900136070872893501182074887119057012393945654193918899145129056911443530288914279128357924753940349347236933442622519850519675741957639564583675552617775213058062276383407118510101599822438563283337500150636725624675182855729196171272117260461324149558239722680824947728130265411739607397255530288824359

echo $N | ./gpu_ecm -save gpuecm.tmp -s $s -n $n $B1
head -n 100 gpuecm.tmp | tail -n 1 > gpuecm2.tmp 

A=`cut -d ";" -f 2 gpuecm2.tmp | cut -c 4-`
echo $N | ecm -save ecm.tmp -A $A -x0 2 $B1 0 > /dev/null
  
t=`cut -d ";" -f 6 ecm.tmp | cut -c 11-`
s=`cut -d ";" -f 6 gpuecm2.tmp |cut -c 11-`

if [ $t -eq $s ] ;then
  echo "Ok! $t"
  rm -f ecm.tmp gpuecm.tmp gpuecm2.tmp
else
  echo "Erreur! $s"
  echo "See ecm.tmp gpuecm2.tmp"
  rm -f gpuecm.tmp
fi

