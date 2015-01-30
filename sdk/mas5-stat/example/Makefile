#
# cvs:affy/sdk/rawq/example/Makefile ---
# $Id: Makefile,v 1.2 2006-01-05 22:01:28 csugne Exp $
#

ROOT = ../..

include ${ROOT}/build/config_vars.mk

SRCS = $(wildcard *.cpp)

OBJS = \
${OBJDIR}/RunAlgorithm.o

INCL += -I../../file
INCL += -I../src
INCL += -I../../rawq

PROGS = ${BINDIR}/mas5-stat-example

include ${ROOT}/build/build_rules.mk

include .depend
