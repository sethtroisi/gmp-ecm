#!/usr/bin/bash

ECM=${1:-./ecm}
DEFAULT_CURVES=$(echo "2^31-1" | ./ecm -gpu 10 | grep -oP '[0-9]*(?= Step 1 took)')
C=${2:-$DEFAULT_CURVES}
B1=128000

if [[ "$C" -ne "$DEFAULT_CURVES" ]]; then
 echo "Changed DEFAULT_CURVES from $DEFAULT_CURVES to $C"
 echo
fi

echo "This script helps find the best \"-gpucurve=<X>\" for your gpu."
echo "$ECM will be run with different multiples of the default ($DEFAULT_CURVES)."
echo "And on four differently sized numbers: 256, 512, 1024, and 2048 bits."
echo "The first line is the CPU timing, then the GPU times for different values of \"-gpucurve\""
echo
echo "Large values tend to produce better throughput but double the time to get the curves."
echo "Using the smallest -gpucurve value within 10% of the best throughput is a good choice."
echo

filtered() {
 echo "$1" | $ECM -v $2 $B1 0 2>&1 | grep -P "CGBN|Step"
}

for number in "(2^269-1)/13822297" "(2^499-1)/20959" "2^997-1" "(2^1877-1)/15017"; do
  echo -e "\n\nTESTING  $number  with B1=$B1"
  filtered "$number"
  echo

  curve_test="$((C / 4)) $((C / 2)) $C $((C * 2)) $((C * 4)) $((C * 8))"
  for curves in $curve_test; do
    filtered "$number" "-v -gpu -gpucurves $curves"
    echo
  done

  B1=$((B1 / 2))
done
