#! /bin/sh

P95_URL=https://mersenne.org/download/software/v30/30.19/p95v3019b12.source.zip
echo "Loading P95 version 30.19 build 12"
[ -d "/tmp/P95/" ] && rm -r /tmp/P95
mkdir /tmp/P95
wget -q -O /tmp/P95/P95_source.zip $P95_URL
unzip -d /tmp/P95 /tmp/P95/P95_source.zip
(cd /tmp/P95/gwnum && make -f make64)
autoreconf -i
./configure --with-gwnum=/tmp/P95/gwnum
make
make check

