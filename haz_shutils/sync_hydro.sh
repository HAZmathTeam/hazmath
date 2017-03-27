#!/bin/sh
# 
# usage is sync
# rsync to another repo which needs hazath libraries and sources

VERSION="0.0.0"

HYDROLENA_HAZ_DIR=../../hydrolena/hazmath-${VERSION}

echo "R-syncing  ${HYDROLENA_HAZ_DIR} on $(date +%Y_%m_%d)..."

set -x
#
## make sure every dir name has backslash at the end of it.
#
rsync -auve --progress --existing --exclude='.git*' \
      --exclude-from='tar.excl' \
      ../ ${HYDROLENA_HAZ_DIR}/

set +x

echo "All synced to ${HYDROLENA_HAZ_DIR}..."
