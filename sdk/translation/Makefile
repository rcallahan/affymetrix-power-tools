#
# cvs:affy/sdk/translation/Makefile ---
#
# $Id: Makefile,v 1.20 2009-04-22 21:13:55 awilli Exp $
#

# sdk_root has a relative path to the top of the tree.
# Once set, suck in our variables and definitions.
sdk_root:=..
include ${sdk_root}/Makefile.defs

# External librarys 
#$(call sdk_set_link_libs,  affysdk affycalvin affyfile affyutil pcrecpp pcre  m)

# These calls define a series of rules used to build executables.
# (files which arent used in the calls will be fed to a default set
# of actions defined in "Makefile2.post".)
#$(call sdk_define_exe,dump-guid,dump-guid.cpp)
$(call sdk_define_exe,apt-dmet-translation,apt-dmet-translation.cpp)
$(call sdk_define_install_exe,apt-dmet-translation)

# While we could use ${sdk_root} to name the file,
# this allows one of the prior makefiles or calls
# to change the "post" makefile by setting the variable.
include ${sdk_makefile_post}


