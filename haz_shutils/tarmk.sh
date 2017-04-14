#!/bin/sh
# 
# usage: ./tarmk.sh 
# Creates a tar file with the latest hazmath sources and examples. 


VERSION="0.0.0"
DNAME="hazmath-${VERSION}"
##FILENAME="hazmath_$(date +%Y_%m_%d)_$VERSION.tar.gz"
###FILENAME="${DNAME}_$(date +%Y_%m_%d).tar.gz"
## for not changing the repo every time, we have a constant name
## independent of the date.

FILENAME="${DNAME}.tar"

FILENAMEz="${FILENAME}.gz"

#set -x

[ -f ${FILENAME} ] && /bin/rm -rf  ${FILENAME}

[ -f ${FILENAMEz} ] && /bin/rm -rf ${FILENAMEz}

echo "Creating TAR file with the current version of HAZMATH: $DNAME on $(date +%Y_%m_%d)"

tar \
    -s:^\.\.:${DNAME}: --exclude-from="tar.excl" \
    -cvf ${FILENAME} ../

#echo "Adding directories needed for output to ${FILENAME}"

#tar --no-recursion \
#    -s:^\.\.:${DNAME}: --exclude-from="tar.excl" \
#    -rvf ${FILENAME} ../examples/*/output

echo "Compressing: ${FILENAME} --> ${FILENAMEz}"

[ -f  ${FILENAME} ]  && gzip -S'.gz' -f ${FILENAME}

#set +x

if [ $? -eq  0 ] ;  then
       echo "*** $FILENAME successfully created. ***"
else
      echo "*** ERROR: ***  the shell exited with a nonzero flag during creation of the archive: ${FILENAME}"
      echo "NO tar file created"
fi
