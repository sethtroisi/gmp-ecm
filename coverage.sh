#!/bin/bash
# script to perform a coverage test
t=$1 # target directory for the html files
d=`mktemp -d /tmp/ecmXXX`
cd $d
svn checkout svn://scm.gforge.inria.fr/svnroot/ecm/trunk ecm
cd ecm
autoreconf -i
./configure --disable-assert
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS00=1
make check VALGRIND=
make clean
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS11=1
make check VALGRIND=
make clean
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS22=1
make check VALGRIND=
make clean
# make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS33=1
# make check VALGRIND=
# make clean
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
make longcheck VALGRIND=
make bench_mulredc
./bench_mulredc
make tune
./tune -v
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm.info ./ --no-external
rm -rf $t
genhtml -o $t/ ecm.info
cd
rm -rf $d
