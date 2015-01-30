#
#

sdk_root:=../../..

include ${sdk_root}/Makefile.defs

# These are four executables 
apg_exe:=apt-snp-model-converter
pgt_exe:=
pgi_exe:=apt-snp-model-converter-quick-test

# This is a variable representing all executables.  
# Used for copying them from output directory for local use.
all_exe:=${apg_exe} ${pgt_exe} ${pgi_exe}

# The variable sdk_exe is a list of the executables we want generated from cpp files in the present directory. 
# The inclusion of the Makefile.defs and Makefile.post ensures that all local .cpp files are
# compiled and that an executable is produced for all values in sdk_exe.  
sdk_exe+=${pgt_exe} ${pgi_exe}

include ${sdk_makefile_post}

copy_executables:${all_exe:%=${sdk_output_bin}/%} 
	${CP} ${all_exe:%=${sdk_output_bin}/%} .

# The following targets are defined in the Makefile.defs file and are invoked when the corresponding
# target without the underscore  '_' is invoked from make.   
_integration_run: copy_executables
	./${pgi_exe}

#_regression_run: copy_executables 
#	./${pgt_exe}

_clean_run:
	-rm -r test-generated


