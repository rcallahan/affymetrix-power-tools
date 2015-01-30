#
# affy/sdk/dbtools/apt-annotation-converter/Makefile ---
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs

# dont put any of these files in the libraries.
$(call sdk_ignore,$(wildcard *.cpp))

# but do compile them all to make this exectuable.
$(call sdk_define_exe,apt-annotation-converter,$(wildcard *.cpp))

#
include ${sdk_makefile_post}
