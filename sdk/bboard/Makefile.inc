#
# /nfs/ss11/harley/Apt/work/apt2-1/affy/sdk/apt2/Makefile.inc ---
#
# $Id: Makefile.inc,v 1.2 2009-11-01 01:03:06 harley Exp $
#

gdb:=/usr/bin/gdb --args
valgrind:=valgrind --tool=memcheck -v --leak-check=full --show-reachable=yes --num-callers=30

# Local Variables:
# mode: makefile
# End:
