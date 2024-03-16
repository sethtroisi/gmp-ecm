#! /bin/sh

P95_URL="https://www.mersenne.org/download/software/v30/30.19/p95v3019b13.source.zip"
echo "Loading P95 version 30.19 build 13"
[ -d "/tmp/P95/" ] && rm -r /tmp/P95
mkdir /tmp/P95
wget -q -O /tmp/P95/P95_source.zip $P95_URL
unzip -d /tmp/P95 /tmp/P95/P95_source.zip
(cd /tmp/P95/gwnum && make -f make64)
autoreconf -i
./configure --with-gwnum=/tmp/P95/gwnum
make
make check
make longcheck
