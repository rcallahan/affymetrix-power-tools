#!/bin/bash
#
# affy/sdk/dox/update-manual-options.sh ---
#
# $Id$
#
# @todo We should add a set of "meta-tags" to our program files
#       which select the files for particuar kinds of processing.
#       A tag to control what is processesed, in additon to the
#       tags used for processing.

# Build a list of cpp files which have "manualOptions",
# but EXCLUDE the chipdesign/ directory as it is not normally built.

# Check we're in the right directory
if [ ! -d ../dox ]; then echo "update-manual-options.sh: should be run from affy/sdk/dox" >2; fi

# set sed command for apple or linux
if [ $(uname) = "Darwin" ]
then
  SED="sed -i \"\""
else
  SED="sed -i"
fi


FILES=$(find .. -name '*.cpp' -or -name '*.dox' | xargs grep -l manualOptions | 
	grep -v -e chipdesign/ | grep -v -e release/ |  sort)

for FILE in $FILES; do

    BASE=${FILE##*/}
    BASE=${BASE%.*}
    APTHELPOUT=${FILE%.*}.help.txt
    APTBIN=../output/`../build/cpucomsys.sh`/bin/$BASE

    export APTBIN
    export APTCPP

    if [ -x "$APTBIN" ]
    then
      echo "Writing $APTHELPOUT"
  
      if [ -e "$APTHELPOUT" ]
      then
        rm $APTHELPOUT
      fi
  
      $APTBIN -help > $APTHELPOUT
  
      if [ $? != 0 ]; then echo "ERROR: Update docs failed on $BASE"; exit 1; fi
  
      ${SED} '/^version:/d' $APTHELPOUT
  
      if [ $? != 0 ]; then echo "ERROR: Update docs failed on $BASE"; exit 1; fi

    else
      echo "ERROR: Update docs could not execute \"$APTBIN -help\""; exit 1;
    fi


done;

# say we are ok.
exit 0
