#!/bin/bash

# Copy required documentation files from dox generated directory to apt/doc
# and sanitize them to be used as online documentation.
#
# TODO: See if this formatting could be done directly by doxygen instead of post-processing.
#

# check called from sdk/dox directory
if [ "$(basename $(pwd -P))" !=  "dox" ]
then
  echo "ERROR: doc-files.sh must be run from the ./sdk/dox directory"
  exit 1
fi

# set manifest file if not set
if [ "${APT_MANIFEST_DOCS}" = "" ]
then
  APT_MANIFEST_DOCS='MANIFEST_docs.txt'
fi

if [ "${APT_REL_PATH}" = "" ]
then
  APT_REL_PATH='../release'
fi

# These are files which go into the binary APT releases
FILES=$(cat ${APT_REL_PATH}/${APT_MANIFEST_DOCS})

DOCTODIR="../apt/doc"
DOCFROMDIR="./html"

# set sed command for apple or linux
if [ $(uname) = "Darwin" ]
then
  SED="sed -i \"\""
else
  SED="sed -i"
fi

# create doc subdirectory if it doesn't exist
if [ ! -d ${DOCTODIR} ]; then mkdir ${DOCTODIR}; fi

for F in $FILES; do
  file_from=${DOCFROMDIR}/$F
  file_to=${DOCTODIR}/$F
  echo "Copying document file '$file_from' to '$file_to'"
  cp  $file_from $file_to
  if [ $? != 0 ]; then 
    echo "failed to find '$F'"
    exit -1
  fi
  if [ ${F##*.} = "html" ]
  then
    echo "Fixing up document file $F"
    # sed substitutions here
    #
    # Remove date from "Generated on ..." line
    # Remove unwanted hrefs to code docs
    # Remove link to pages.html
    # replace doxygen header tabs with APT tabs
    echo "Fixing tabs in file $F"
    ${SED} \
    -e 's/Generated on .* for Affymetrix Power Tools/Generated for Affymetrix Power Tools/g' \
    -e 's/<a class=\"el\" href=\"class.*>\(.*\)<\/a>/\1/g' \
    -e 's/<a .*href=\".*_8h\.html\".*>\(.*\)<\/a>/\1/g' \
    -e 's/<a .*href=\".*_8cpp.*\.html\".*>\(.*\)<\/a>/\1/g' \
    -e '/<a .*href=\"pages\.html\".*>/ d' \
    -e '
       /<div [^>]*class=\"tabs\">/,/<\/div>/ {
	  /^[ \t]*<li[ >]/ d
	  /<ul class=\"tablist\">/a\
             <li><a href=\"index.html\"><span>Main&nbsp;Page</span></a></li>\
             <li><a href=\"CHANGE.html\"><span>Change&nbsp;Log</span></a></li>\
             <li><a href=\"VIGNETTE.html\"><span>Vignettes</span></a></li>\
             <li><a href=\"FAQ.html\"><span>FAQs</span></a></li>\
             <li><a href=\"LICENSE.html\"><span>License</span></a></li>\
             <li><a href=\"FILE-FORMATS.html\"><span>File Formats</span></a></li>\
             <li><a href=\"PLATFORMS.html\"><span>Platforms</span></a></li>
	}
    ' \
    -e '/^<div id=\"titlearea\">/,/^<\/div>/ d' \
    -e '/<div class=\"title\">\([^<]*\)<\/div>/s//<h1>\1<\/h1>/' \
    ${DOCTODIR}/$F
    if [ $? != 0 ]; then 
          echo "Failed to fix up file $F"
          exit -1
    fi
    # If file is linked from header tabs then make status of corresponding tab current
    ${SED} -e '/<li><a href=\"'"$F"'\">/ s/<li>/<li class=\"current\">/' ${DOCTODIR}/$F
    if [ $? != 0 ]; then 
          echo "Failed to fix up file $F"
          exit -1
    fi
    rm ${DOCTODIR}/${F}.bak
  fi
done;

# Search copied HTML and CSS files for included .png files and copy them over too.
#
PNGS="$(cat ${DOCTODIR}/*.html | sed 's/<img /\n<img /g' | sed -n 's/.*src=\"\([^ ]*\.png\)\".*/\1/p' | sort | uniq)"
PNGS="$(cat ${DOCTODIR}/*.css | sed -n 's/.*[^_a-zA-Z0-9]\([_a-zA-Z0-9]*\.png\).*/\1/p' | sort | uniq) $PNGS"
for F in $PNGS; do
    echo "Copying document file '$F'"
    cp ${DOCFROMDIR}/$F ${DOCTODIR}/$F
    if [ $? != 0 ]; then 
        echo "Failed to find '$F'"
    fi
done;

echo "Copying over AGCC/GCOS docs"
if [ ! -d ${DOCTODIR}/gcos-agcc ]; then mkdir ${DOCTODIR}/gcos-agcc; fi;
cp -r ${DOCFROMDIR}/gcos-agcc/* ${DOCTODIR}/gcos-agcc/

exit 0
