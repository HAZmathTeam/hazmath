#!/bin/sh
# 
# usage: ./tarmk.sh 
# Creates a tar file with the latest hazmath sources and examples. 


VERSION="1.0.0"
DNAME="hazmath-${VERSION}"
##FILENAME="hazmath_$(date +%Y_%m_%d)_$VERSION.tar.gz"
###FILENAME="${DNAME}_$(date +%Y_%m_%d).tar.gz"
## for not changing the repo every time, we have a constant name
## independent of the date.

FILENAME="${DNAME}.tar"

FILENAMEz="${FILENAME}.gz"

#set -x

[ -f ../tarball/${FILENAME} ] && /bin/rm -rf  ../tarball/${FILENAME}

[ -f ../tarball/${FILENAMEz} ] && /bin/rm -rf ../tarball/${FILENAMEz}

if [ -d ../tarball/${DNAME} ] ; then
    echo "${DNAME}/ exists"
else
    mkdir ../tarball/${DNAME}
fi
    rsync -avv \
	  --exclude-from='tar.excl' \
	  ../ ../tarball/${DNAME}/
rsync -avv \
      --include-from='grids.incl' \
      --exclude='*' \
      ../ ../tarball/${DNAME}/


echo "Creating TAR file with the current version of HAZMATH: $DNAME on $(date +%Y_%m_%d)"

cd ../tarball ; tar  cvf ${FILENAME} ${DNAME} && /bin/rm -rf ${DNAME}

#set +x
###tar -s:^\.\.:${DNAME}: --exclude-from="tar.excl" -cvf ${FILENAME} ../

## create a symbolic link.

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
