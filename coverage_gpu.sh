#!/bin/bash
# script to perform a coverage test, with gpu support

t=$1 # target directory for the html files
d=`mktemp -d /tmp/ecmXXX`
cd $d
git clone https://gitlab.inria.fr/zimmerma/ecm.git
cd ecm
autoreconf -i
# ./configure --enable-gpu=sm_30 --with-gmp=/users/caramel/logiciels/gmp-6.0.0/core2/ --with-cuda=/usr/local/cuda-5.0.old/ --with-cc-for-cuda=/users/caramel/logiciels/gcc-4.3.6/x86_64/bin/
# ./configure --enable-gpu=sm_30 --with-cuda=/tmp/cuda
./configure --enable-gpu --disable-shared --enable-static
# make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
make
export LD_LIBRARY_PATH=/usr/lib/cuda/lib64:.
./test.gpuecm ./ecm
make longcheck VALGRIND=
make bench_mulredc
./bench_mulredc
make tune
./tune -v
# geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm.info ./ --no-external
# rm -rf $t
# genhtml -o $t/ ecm.info
cd
rm -rf $d
