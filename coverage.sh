#!/bin/bash
# script to perform a coverage test
t=$1 # target directory for the html files
d=`mktemp -d /tmp/ecmXXX`
cd $d
svn checkout svn://scm.gforge.inria.fr/svnroot/ecm/trunk ecm
cd ecm
autoreconf -i
./configure
make CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
make longcheck
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename ecm.info ./ --no-external
rm -rf $t
genhtml -o $t/ ecm.info
cd
rm -rf $d
