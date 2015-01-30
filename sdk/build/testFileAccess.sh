#!/bin/bash
#
# affy/sdk/testFileAccess.sh ---
#
# $Id$
#
# Used to identify files accessed by regressionn tests
#

VERB=0

# Get Options
while getopts "d:o:t:r:v?" opt; do
    case $opt in
        "d")
            TESTDIR=$OPTARG
            ;;
        "o")
            OUTDIR=$OPTARG
            ;;
        "t")
            TEST=$OPTARG
            ;;
        "r")
            RFDIR=$OPTARG
            ;;
        "v")
	    VERB=1
	    ;;
        "?")
            echo
            echo "usage: testFileAccess.sh"
            echo
	    echo "   -d test directory (in which 'make -k regression' will be run)"
	    echo "   -o output directory (default .)"
	    echo "   -t test name (default test dir name)"
            echo "   -r regression files path (default ./regression-data/data/cel)"
            exit 0
    esac
done

# Set up variables and check files
#

if [ "$TESTDIR" = "" ]
then
  TESTDIR=$(pwd -P)
fi
if [ $VERB == 1 ]; then echo "Setting test directory to $TESTDIR"; fi

if [ "$OUTDIR" = "" ] 
then
  OUTDIR=$(pwd -P)
fi
if [ $VERB == 1 ]; then echo "Setting output directory to $OUTDIR"; fi

if [ "$TEST" = "" ] 
then
  TEST=${TESTDIR##*/}
fi
if [ $VERB == 1 ]; then echo "Setting test name to $TEST"; fi

if [ "$RFDIR" = "" ] 
then
  RFDIR=./regression-data/data/cel
fi
if [ $VERB == 1 ]; then echo "Setting regression data directory to $RFDIR"; fi

# set time stamp

TSFILE="$OUTDIR/$TEST.ts"

touch "$TSFILE"

# run make

(cd "$TESTDIR"; make -k regression)

# check for accessed files

if [ $VERB == 1 ]; then echo "(cd \"${RFDIR}\"; find . -type f -anewer \"${TSFILE}\" | sort > \"${OUTDIR}/${TEST}_files\")"; fi
(cd "$RFDIR"; find . -type f -anewer "$TSFILE" | sort > "${OUTDIR}/${TEST}_files")

#rm "$TSFILE"


