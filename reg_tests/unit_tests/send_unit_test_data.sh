#!/bin/bash

# Variables
INOFILE=.info
VERSIONFILE=.version

# This script upload data for JOREK unit tests.
UNITTESTDIR="$1"
if [ -z "${DAV_URL}" ]; then
    DAV_URL="http://jorek.eu/dav_nrt/"
fi
if [ -z "${UNITTESTDIR}" ]; then
    printf "Failed: requires directory name!\n"
    exit 1 
fi

cd  ${UNITTESTDIR}
# obtain path to the unit test directory
unittestcasedir=`readlink -f ${PWD}`
# extract all files
mapfile -t unit_test_files <<< "$(ls -1)"

# check if folder is empty
if [ -z "${unit_test_files[@]}" ]; then
    printf "Failed: empty folder!\n"
    exit 1
fi

#compute the checksum for all files
VERSION="_`cat ${unit_test_files[@]} | md5sum | sed -e 's/^ //' -e 's/ .*$//'`"
echo ${VERSION} > ${VERSIONFILE}

printf "date: "      > ${INFOFILE}
date                >> ${INFOFILE}
printf "hostname: " >> ${INFOFILE}
hostname -f         >> ${INFOFILE}
printf "branch: "   >> ${INFOFILE}
git rev-parse --abbrev-ref HEAD >> ${INFOFILE}
printf "user: "     >> ${INFOFILE}
echo $USER          >> ${INFOFILE}

# generating the name of the unit test tarball
unittestname=$(basename $unittestcasedir)
TGZFILE=${unittestname}${VERSION}.tgz

echo "Creating tarball ${TGZFILE}"
tar cvzf ${TGZFILE} ${INFOFILE} ${unit_test_files[@]} || exit 1
 
echo "Uploading ${TGZFILE}"
UNITTESTNAME=$(basename $UNITTESTDIR)
curl -s -u nrt:nrt_21745XtL -T ${TGZFILE} ${DAV_URL}
if [ $? -eq 0 ]; then
  printf "Success\n"
  exit 0
else
  printf "Failed\n"
  exit 1
fi
