#!/usr/bin/bash

ECM=${1:-./ecm}
DEFAULT_CURVES=$(echo "2^31-1" | ./ecm -gpu 10 | grep -oP '[0-9]*(?= parallel curves)')
C=${2:-$DEFAULT_CURVES}
B1=64000

if [[ "$C" -ne "$DEFAULT_CURVES" ]]; then
 echo "Changed DEFAULT_CURVES from $DEFAULT_CURVES to $C"
fi

filtered() {
 echo "$1" | $ECM -v $2 $B1 0 2>&1 | grep -P "Input|CGBN|Step"
}

for number in "(2^269-1)/13822297" "(2^499-1)/20959" "2^997-1" "2^1907-1"; do
  echo -e "\n\nTESTING $number B1=$B1"
  filtered "$number"
  filtered "$number" "-gpu"
  echo

  curve_test="$((C / 4)) $((C / 2)) $C $((C * 2)) $((C * 4))"
  for curves in $curve_test; do
    filtered "$number" "-cgbn -gpucurves $curves"
    echo
  done

  B1=$((B1 / 2))
done
