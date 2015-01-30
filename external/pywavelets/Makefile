#
# svn://svn.ev.affymetrix.com/apt/externals/pywavelets/trunk/Makefile
#

# this should be set already, but make sure it is.
ifndef CPUCOMSYS
export CPUCOMSYS:=$(shell ../../sdk/build/cpucomsys.sh)
endif

#
default: all

#
AR=ar rc
RANLIB=ranlib

#
OBJS:=convolution.o wavelets.o wt.o common.o
DOT_H:=common.h convolution.h wavelets.h wavelets_coeffs.h

.PHONY: all install apt-install clean

all: libpywavelets.a install

libpywavelets.a: ${OBJS}
	$(AR) ${@} $(OBJS)
	$(RANLIB) $@

%.o: %.c ${DOT_H}
	${CXX} -fPIC ${CFLAGS} -c ${<} -o ${@}

clean:
	-rm ${OBJS} libpywavelets.a

install apt-install: libpywavelets.a
	cp libpywavelets.a ../../sdk/output/${CPUCOMSYS}/lib
