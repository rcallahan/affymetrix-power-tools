# This Makefile is for use within the APT package
#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.18 2009-09-10 21:16:01 awilli Exp $
#

#
sdk_root:=..
include ${sdk_root}/Makefile.defs
# See the note in this file.
include ./Makefile.inc

#
sdk_subdirs:=\
	exception/src parameter/src utils/src template/src parsers/src \
	writers/src fusion/src data/src \
	exception/test parameter/test utils/test parsers/test \
	writers/test fusion/test data/test \
    portability/test \
	Extractor/src converters array/src array/test
#
include ${sdk_makefile_post}
