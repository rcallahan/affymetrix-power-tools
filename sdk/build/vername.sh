#!/bin/bash

### vername.sh - generate svn version or release version name
#
# For SVN revision:
#   First check if VERSION file contains SVNREV line.
#   If that doesn't work then use svnversion command.
# For SVN URL:
#   First check if VERSION file contains SVNURL line.
#   If that doesn't work then use svn info
# For release version:
#   Check VERSION file for RELEASE line
#   Otherwise use entire contents of VERSION file
#

# Parse parameters
if [ "$#" -eq 0 ]
then
  echo "Usage: vername.sh (-r|-s) [-v]" 1>&2
  exit 1
fi

while [ "${1}" != "" ]
do
  if [ "${1}" = "-v" ]
  then
    opt_verbose=1
  elif [ "${1}" = "-s" ]
  then
    opt_svn=1
  elif [ "${1}" = "-u" ]
  then
    opt_url=1
  elif [ "${1}" = "-r" ]
  then
    opt_release=1
  else
    echo "Usage: vername.sh (-r|-s) [-v]" 1>&2
    exit 1
  fi
  shift
done

# should be called from sdk or a subdirectory
if [ "${sdk_root}" = "" ]
then
  # look up the tree for a version file.
  dir=$(pwd)
  while [ ! -f ${dir}/VERSION ]
  do
    dir=${dir%/*}
    if [ "${dir}" = "" ]
    then
      echo "Unable to find 'VERSION' above ." 1>&2
      exit 1
    fi        
    if [ "$opt_verbose" = 1 ]
    then
      echo "${0}: search ${dir}/VERSION" 1>&2
    fi        
  done
  version_file=${dir}/VERSION
else
  # sdk_root is set in env, so use it directly.
  version_file=${sdk_root}/../VERSION
fi

#
if [ ! -e "${version_file}" ]
then
  echo "Unable to find 'VERSION' above 'sdk'." 1>&2
  exit 1
fi

if [ "$opt_svn" = 1 ]
then
  svn_ver=$(awk '/^SVNREV/ {print $2}' ${version_file})
  if [ "$svn_ver" = "" ]
  then
    # check for svnversion
    which svnversion > /dev/null 2>&1
    if [ $? = 0 ]; then
      svn_ver=$(svnversion)
    else
      svn_ver='exported'
    fi
  fi
  echo "$svn_ver"

elif [ "$opt_url" = 1 ]
then
  svn_url=$(awk '/^SVNURL/ {print $2}' ${version_file})
  if [ "$svn_url" = "" ]
  then
    # check for svn
    which svn > /dev/null 2>&1
    if [ $? = 0 ]; then
        svn_url=$(svn info "${version_file}" | grep 'URL' | sed -e 's/URL: //' -e 's/\/VERSION//')
    else
        svn_url=""
    fi
  fi
  echo "$svn_url"

else # release version

  ver_name=$(awk '/^RELEASE/ {print $2}' ${version_file})
  if [ "$ver_name" = "" ]
  then
    ver_name=$(head -n 1 ${version_file})
  fi
  echo "$ver_name"
fi

exit 0
