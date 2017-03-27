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

[ -f ${FILENAME}  ] && /bin/rm -rf ${FILENAME}

echo "Creating TAR file with the current version of HAZMATH: $DNAME on $(date +%Y_%m_%d)"

#set -x

tar \
    -s:^\.\.\/:${DNAME}/: --exclude-from="tar.excl" \
    -zcvf ${FILENAME} ../

#set +x

if [ $? -eq  0 ] ;  then
       echo "*** $FILENAME successfully created. ***"
else
      echo "*** ERROR: ***  the shell exited with a nonzero flag during creation of the archive: ${FILENAME}"
      echo "NO tar file created"
fi
