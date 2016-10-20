#!/bin/bash
# script to perform a coverage test
#
# In case of problems, make sure gcov has the same version number as gcc.
# If not, you might want to add the flag
#     --gcov-tool <name_of_the_correct_gcov>
# in geninfo. 
# Also you may add CC=<name_of_the_correct_gcc> to the make command.
#
t=$1 # target directory for the html files
d=`mktemp -d /tmp/ecmXXX`
cd $d
svn checkout svn://scm.gforge.inria.fr/svnroot/ecm/trunk ecm
cd ecm
autoreconf -i
./configure --disable-assert
echo "Testing PARAMS00"
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS00=1
./test.ecm ./ecm
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm00.info ./ --no-external
make clean
echo "Testing PARAMS11"
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS11=1
./test.ecm ./ecm
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm11.info ./ --no-external
make clean
echo "Testing PARAMS22"
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS22=1
./test.ecm ./ecm
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm22.info ./ --no-external
make clean
echo "Testing PARAMS33"
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage" PARAMS33=1
./test.ecm ./ecm
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm33.info ./ --no-external
make clean
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
echo "Testing longcheck"
make longcheck VALGRIND=
echo "Testing bench_mulredc"
./bench_mulredc -v
echo "Testing tune"
./tune -v
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm.info ./ --no-external
rm -rf $t
genhtml -o $t/ ecm.info ecm00.info ecm11.info ecm22.info ecm33.info
cd
rm -rf $d
