#
# cvs: ---
#
# $Id: Makefile,v 1.11 2009-09-30 18:39:53 csugne Exp $
#

#

rm:=/bin/rm -rf
out_dir:=output
sdk_root:=../..

atd_exe:=apt-dmet-translation
reg_exe:=dmet-translation-test

include ${sdk_root}/Makefile.defs
#
sdk_cpp_ignore=$(wildcard *.cpp)
_integration_run+=_mk_output_dir
_integration_run+=integration_1 
_clean_run+=clean_integration
# Add the output directory to the list.
_mk_output_dir: | output
#

# External librarys 
#$(call sdk_set_link_libs, affysdk affycalvin affyfile affyutil pcrecpp pcre  m)

$(call sdk_define_exe,${reg_exe},$(wildcard *.cpp))

include ${sdk_makefile_post}

output:
	mkdir ${out_dir}

_integration_exe: ${sdk_output_bin}/${atd_exe}

cp_pgms: ${sdk_output_bin}/${reg_exe} ${sdk_output_bin}/${atd_exe}
	${CP} ${sdk_output_bin}/${atd_exe} ${sdk_output_bin}/${reg_exe} .

integration_1: clean_integration cp_pgms
	mkdir -p ${out_dir}
	./${reg_exe}

integration_test1: 
	./${reg_exe} -i input1

integration_test1_debug:
	./${atd_exe} -i input1

integration_test2: 
	./${reg_exe} -i input2

integration_test2_debug: 
	./${atd_exe} -i input2

integration_test3: 
	./${reg_exe} -i input3

integration_test3_debug: 
	./${atd_exe} -i input3

clean_integration:
	-${rm} ${out_dir}
	mkdir -p ${out_dir}
