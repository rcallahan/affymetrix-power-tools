#!/bin/sh

### figure out what kind of machine we are on.

uname=`uname 2>/dev/null`
uname_o=`uname -o 2>/dev/null`

if [ "$uname" = "Darwin" ]; then
	uname_p=`uname -p 2>/dev/null`
	uname_m=`uname -m 2>/dev/null`
        uname_r=`uname -r 2>/dev/null`
	uname_r_major=${uname_r%%.*}
	if [ "$uname_r_major" = "11" ]; then
	        echo "${uname_m}-apple-lion"
	elif [ "$uname_r_major" = "10" ]; then
		echo "${uname_p}-apple-snowleopard"
	elif [ "$uname_r_major" = "9" ]; then
		echo "${uname_p}-apple-leopard"
	elif [ "$uname_r_major" = "8" ]; then
		echo "${uname_p}-apple-tiger"
	elif [ "$uname_r_major" = "7" ]; then
		echo "${uname_p}-apple-panther"
	else
		echo "unknown" # -darwin-system
	fi
elif [ "$uname" = "SunOS" ]; then
	uname_p=`uname -p 2>/dev/null`
	echo "${uname_p}-sun-solaris"
elif [ "$uname_o" = "Cygwin" ]; then
	echo "i386-intel-cygwin"
elif [ "$uname" = "Linux" ]; then
        if [ "$(uname -v | grep 'Ubuntu')" != "" ]
        then
          linuxtype="ubuntu"
        else
          linuxtype="linux"
        fi
        # uname -i doesnt work on ubuntu, returns i386 on Fedora
        # uname -m returns i686 on Fedora
	uname_m=`uname -m 2>/dev/null`
	if [ "$uname_m" = "i386" ] || [ "$uname_m" = "i686" ] ; then
		echo "i386-intel-${linuxtype}"
	elif [ "$uname_m" = "x86_64" ]; then
		echo "x86_64-intel-${linuxtype}"
	elif [ "$uname_m" = "ppc" ]; then
		echo "powerpc-apple-${linuxtype}"
	else
		echo "unknown"
	fi
else
	echo "unknown"
fi

# There is also "i386-intel-win32" which is not covered here.


