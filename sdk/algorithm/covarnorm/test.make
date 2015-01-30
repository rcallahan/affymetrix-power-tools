sdk_root:=../../..
include ${sdk_root}/Makefile.defs
#                  
sdk_cpp_cflags+= -I.. ${newmat_cflags}
sdk_cpp_lflags+=${newmat_lflags}
#                  
$(call sdk_define_exe,toy.norm,toy.norm.cpp)
#                  
include ${sdk_makefile_post}
