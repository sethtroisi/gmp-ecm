#! /bin/sh

P95_URL=https://www.mersenne.org/download/software/v30/30.8/p95v308b15.source.zip
[ -d "/tmp/P95/" ] && rm -r /tmp/P95
mkdir /tmp/P95
wget -q -O /tmp/P95/P95_source.zip $P95_URL
unzip -d /tmp/P95 /tmp/P95/P95_source.zip
(cd /tmp/P95/gwnum && make -f make64)
autoreconf -i
./configure --with-gwnum=/tmp/P95/gwnum
make
make check

