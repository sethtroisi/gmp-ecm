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

N=780851074980451083207022065351299001484424461908031846796515069313846989720932551324409454457573724028694607215099661331321033861266619279756898656992234612627344143511823522694189489558298076363877298043469771509953319662285499302732028963807960560858768120637464833635395868634144783500288962345010601

rm -f gpuecm.tmp gpuecm2.tmp ecm.tmp

echo $N | ./gpu_ecm -save gpuecm.tmp -s $s -n $n $B1
index=`expr $RANDOM % 448` 
head -n $index gpuecm.tmp | tail -n 1 > gpuecm2.tmp 


A=`cut -d ";" -f 2 gpuecm2.tmp | cut -c 4-`
echo $N | ecm -save ecm.tmp -A $A -x0 2 $B1 0 > /dev/null
  
t=`cut -d ";" -f 6 ecm.tmp | cut -c 11-`
s=`cut -d ";" -f 6 gpuecm2.tmp |cut -c 11-`

if [ $t -eq $s ] ;then
  echo "Ok! $t"
  rm -f ecm.tmp gpuecm.tmp gpuecm2.tmp
else
  echo "Error with lines $index (find $s, expected $t)"
  echo "See ecm.tmp gpuecm2.tmp"
  rm -f gpuecm.tmp
fi

