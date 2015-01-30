#!/bin/bash

(
    find .. -name \*.dox; 
    find .. -name \*.h; 
    find .. -name \*.cpp
) | sort | ./doc-find-section-conflicts -stdin > /var/tmp/$$.apt-section-conflicts

WC=`wc -l /var/tmp/$$.apt-section-conflicts | awk '{print $1}'`
if [ $WC -ne 0 ]; then
    cat /var/tmp/$$.apt-section-conflicts
    rm /var/tmp/$$.apt-section-conflicts
    exit -1
else
    exit 0
fi
