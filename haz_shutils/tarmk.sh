#!/bin/sh
# 
# usage is tarmk.sh 
# Edit to include various directories and files into tarball that will be exhanged to users

VERSION="0.0.0"
DNAME="hazmath-${VERSION}"
##FILENAME="hazmath_$(date +%Y_%m_%d)_$VERSION.tar.gz"
###FILENAME="${DNAME}_$(date +%Y_%m_%d).tar.gz"
## for not changing the repo every time, we have a constant name
## independent of the date.
FILENAME="${DNAME}.tar.gz"

echo "Creating Files for $DNAME on $(date +%Y_%m_%d)"

set -x

tar \
    -s:^\.\.\/:${DNAME}/: --exclude-from="tar.excl" \
    -zcvf ${FILENAME} ../

set +x

echo "$FILENAME created."
