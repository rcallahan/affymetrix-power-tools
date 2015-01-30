#
# affy/sdk/bboard/bboard/Makefile ---
#
# $Id: Makefile,v 1.3 2009-11-01 01:03:07 harley Exp $
#

#
sdk_root:=../..
include ${sdk_root}/Makefile.defs
include ${sdk_root}/bboard/Makefile.inc
#
sdk_cpp_ignore:=Apt2Main.cpp apt2-main-test1.cpp apt2-main-test2.cpp 
#
$(call sdk_define_exe,apt-run-pnode,apt-run-pnode.cpp)
#
include ${sdk_makefile_post}

#
apt_run_pnode_exe:=${sdk_output_bin}/apt-run-pnode

_test: _test_arp

_test_arp: _test_arp_null _test_arp_vecstat

_test_arp_null: ${apt_run_pnode_exe}
	${apt_run_pnode_exe} inputs/arp-1-i.bb PN_Null arp-1-o.bb

_test_arp_vecstat: ${apt_run_pnode_exe}
	-rm arp-{2,3}-o.bb
	${apt_run_pnode_exe} inputs/arp-2-i.bb PN_TestVecStat arp-2-o.bb
	${apt_run_pnode_exe} inputs/arp-3-i.bb PN_TestVecStat arp-3-o.bb
	cat arp-{2,3}-o.bb

#
APT2_MAIN_TEST1_CPP:= apt2-main-test1.cpp Apt2Main.cpp
./apt2-main-test1: ${APT2_MAIN_TEST1_CPP}
	g++ -o ${@} ${APT2_MAIN_TEST1_CPP}

_test_apt2_main_test1: ./apt2-main-test1
# a simple test.
	./apt2-main-test1 arg1 arg2 arg3
# this should error out.
	./apt2-main-test1
