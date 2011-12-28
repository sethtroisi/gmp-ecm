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

./gpu_ecm -save gpuecm.tmp -s $s -n $n $B1 < c306
head -n 100 gpuecm.tmp | tail -n 1 > gpuecm2.tmp 

A=`cut -d ";" -f 2 gpuecm2.tmp | cut -c 4-`
ecm -save ecm.tmp -A $A -x0 2 $B1 0 < c306 > /dev/null
  
t=`cut -d ";" -f 6 ecm.tmp | cut -c 11-`
s=`cut -d ";" -f 6 gpuecm2.tmp |cut -c 11-`

if [ $t -eq $s ] ;then
  echo "Ok! $t"
else
  echo "Erreur! $s"
fi

rm -f ecm.tmp gpuecm.tmp gpuecm2.tmp
