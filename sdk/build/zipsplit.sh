#!/bin/bash
#
# zipsplit.sh
#
# Make a series on <1G zip files from specified files
# 
# Usage: zipsplit.sh <output-file-prefix> <files>
#

# First check sufficient arguments
if [[ $1 == '' || $2 == '' ]]
then
  echo "Usage: zipsplit.sh <output-file-prefix> <files>"
  exit 1
fi

# set constants
VERB=1
MAXSIZE=$((1024 * 1024 * 1024))


# read file prefix
outpref=$1
shift
if [ VERB != 0 ]; 
then 
  echo "OUTPUT PREFIX = $outpref"
  echo "MAXSIZE = $MAXSIZE"
fi

# set up variables
zipnum=1;                               # number of current zip file
zipfile="${outpref}.${zipnum}.zip";     # current zip file name
filelist="";                            # current list of files to be zipped
avail=${MAXSIZE};                       # available space for current zip file
total=0;                                # total size of current files in $filelist

# Read remaining arguments and process

while [[ $1 != '' ]]
do
  # read file and check it's valid
  file=$1
  if [[ ! -f "${file}" && ! -d "${file}" ]]
  then
    echo "FILE ${file} IS NOT VALID" >&2
    exit 1
  fi
  # get size
  size=$(du -b "${file}" | cut -f 1)
  if [ VERB != 0 ]; 
  then 
    echo "FILE: ${file} SIZE: ${size}   TOTAL: ${total}"
  fi
  # check total with curretn file greater than available space
  if [[ $(( ${total} + ${size} )) -gt ${avail} ]]
  then
    # zip files
    if [ VERB != 0 ]; 
    then 
      echo "ADDING FILES ${filelist} TO ${zipfile}"
      echo "TOTAL SIZE OF FILE LIST = ${total}"
      echo "AVAILABLE SIZE =          ${avail}"
    fi
    zip "$zipfile" ${filelist}
    #check size of zip file
    zfsize=$(du -b "$zipfile" | cut -f 1)
    avail=$((${MAXSIZE} - ${zfsize}))
    if [ VERB != 0 ]; 
    then 
      echo "NEW SIZE OF ${zipfile} = ${zfsize}"
      echo "NEW AVAILABLE SIZE = ${avail}"
      echo "SIZE OF NEXT FILE =  ${size}"
    fi
    #
    # start new zip file if necessary
    # make a temp zip of $1 to see if it's small enough to fit in current zipfile
    tmpzip="/tmp/tz.$$.zip"
    rm "$tmpzip"
    zip "$tmpzip" "${file}"
    tzsize=$(du -b "$tmpzip" | cut -f 1)
    if [ VERB != 0 ]; 
    then 
      echo "SIZE OF NEXT ZIP =   ${tzsize}"
    fi
    if [[ ${avail} -lt ${tzsize} ]]
    then
      zipnum=$((zipnum + 1))
      zipfile="${outpref}.${zipnum}.zip"; 
      avail=${MAXSIZE};
      if [ VERB != 0 ]; 
      then 
        echo "STARTING NEW ZIPFILE: ${zipfile}"
      fi
    fi
    # reset file list
    filelist=""
    total=0
  fi
  # add to file list
  filelist="${filelist} ${file}"
  total=$((${total} + ${size}));

  # next argument
  shift
done

# zip remaining files
if [ VERB != 0 ]; 
then 
  echo "ADDING FILES ${filelist} TO ${zipfile}"
  echo "TOTAL SIZE OF FILE LIST = ${total}"
  echo "AVAILABLE SIZE =          ${avail}"
fi
zip "$zipfile" ${filelist}

# make list of files
(for f in $(ls "${outpref}".*.zip); do unzip -l "$f"; done) | awk '{print $4}' | grep -v -e '^Name$' -e '^-*$' | sort > "${outpref}.files"
