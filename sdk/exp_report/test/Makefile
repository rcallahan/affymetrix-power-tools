#
# cvs:affy/sdk/file/Makefile ---
#
# $Id: Makefile,v 1.3 2009-09-25 17:48:59 mspald Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
#
sdk_cpp_ignore=$(wildcard *.cpp)
#
#ExpressionReportWorkflowTest.cpp
testfiles=\
    ChangeStatsTest.cpp \
    DetectionStatsTest.cpp \
    ExpressionControlResultTest.cpp \
    ExpressionControlsParameterExtractionTest.cpp \
    ExpressionProbeSetFileExtractionTest.cpp \
    ExpressionProbeSetReporterTest.cpp \
    ExpressionRPTFileDataTest.cpp \
    ExpressionReportControlsTest.cpp \
    ExpressionReportDataTest.cpp \
    ProbeSetStatsTest.cpp \
    ReportTest.cpp \
#
$(call sdk_define_check_exe,test-exp-report,$(testfiles))
#
include ${sdk_makefile_post}
