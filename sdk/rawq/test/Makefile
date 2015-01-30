#
# cvs:affy/sdk/chipstream/CPPTest/Makefile
# $Id: Makefile,v 1.3 2009-09-25 17:49:06 mspald Exp $
#

ROOT = ../..
include ${ROOT}/build/config_vars.mk

# Source files for dependencies.
SRCS = $(wildcard *.cpp)


RAWQTEST = ${TEST_BINDIR}/RawQTest
RAWQWORKFLOWTEST = ${TEST_BINDIR}/RawQWorkflowTest

TEST_OBJS = \
	${TEST_OBJDIR}/RawQTest.o \
	${TEST_OBJDIR}/RawQWorkflowTest.o \
	${TEST_OBJDIR}/CPPTest.o 


TEST_PROGS = ${SELFCREATETEST} \
	${RAWQTEST} \
	${RAWQWORKFLOWTEST}

include ${ROOT}/build/build_rules.mk

check:: ${TEST_OBJS} libtest

libtest:  ${TEST_PROGS}
	${SELFCREATETEST}
	${RAWQTEST}
	${RAWQWORKFLOWTEST}

include .depend


