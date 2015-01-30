#!/bin/bash

### generate name for distribution version
#
# If VERSION has specific version number, then "apt-<VERSION>"
# otherwise "apt-anon-<DATE>"

if [ "${1}" = "-v" ]
then
  shift
  opt_verbose=1
fi

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
      echo "${0}: search ${dir}/VERBOSE" 1>&2
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

ver_name=$(awk '/^RELEASE/ {print $2}' ${version_file})
if [ "$ver_name" = "" ]
then
  ver_name=$(head -n 1 ${version_file})
fi

if [ "$ver_name" = "NON-OFFICIAL-RELEASE" ]
then
  echo "apt-anon-$(date +%Y%m%d)"
else
  echo "apt-${ver_name}"
fi
