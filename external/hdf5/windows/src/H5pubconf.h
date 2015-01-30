
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* H5pubconf.h is adapted from UNIX platform and manually maintained on the windows platform. */

/*#define H5_HAVE_TM_ZONE 1 windows do not use this constant.*/  
#define H5_MALLOC_WORKS 1

/* code warrior returns 0 in malloc(0) */
#if defined(__MWERKS__)
#undef H5_MALLOC_WORKS 
#endif

/* 
code warrior v.8 does not allow shared writing by default; 
the feature can be enabled by defining
_MSL_ALLOW_SHARED_WRITING to 1
in the file file_io.win32.c and including it on the projects
*/
#if defined(__MWERKS__)
#define H5_NO_SHARED_WRITING
#endif


#define H5_STDC_HEADERS 1
#define H5_HAVE_ATTRIBUTE 1
#undef H5_HAVE_ATTRIBUTE
#define H5_HAVE_LARGE_HSIZET 1 
#ifdef __MWERKS__ 
#define H5_PRINTF_LL_WIDTH "ll" 
#else 
#define H5_PRINTF_LL_WIDTH "I64" 
#endif
#define H5_HAVE___int64
#define H5_SIZEOF___INT64 8 
#define H5_SIZEOF_CHAR 1
#define H5_SIZEOF_DOUBLE 8
#define H5_SIZEOF_FLOAT 4
#define H5_SIZEOF_INT 4
#define H5_SIZEOF_INT16_T 0
#define H5_SIZEOF_INT32_T 0
#define H5_SIZEOF_INT64_T 0
#define H5_SIZEOF_INT8_T 0
#define H5_SIZEOF_INT_FAST16_T 0
#define H5_SIZEOF_INT_FAST32_T 0
#define H5_SIZEOF_INT_FAST64_T 0
#define H5_SIZEOF_INT_FAST8_T 0
#define H5_SIZEOF_INT_LEAST16_T 0
#define H5_SIZEOF_INT_LEAST32_T 0
#define H5_SIZEOF_INT_LEAST64_T 0
#define H5_SIZEOF_INT_LEAST8_T 0
#define H5_SIZEOF_LONG 4
#define H5_SIZEOF_LONG_LONG 0
#define H5_SIZEOF_LONG_DOUBLE 8
#define H5_SIZEOF_OFF_T 4
#define H5_SIZEOF_SHORT 2
#define H5_SIZEOF_SIZE_T 4
#define H5_SIZEOF_SSIZE_T 0
#define H5_SIZEOF_UINT16_T 0
#define H5_SIZEOF_UINT32_T 0
#define H5_SIZEOF_UINT64_T 0
#define H5_SIZEOF_UINT8_T 0
#define H5_SIZEOF_UINT_FAST16_T 0
#define H5_SIZEOF_UINT_FAST32_T 0
#define H5_SIZEOF_UINT_FAST64_T 0
#define H5_SIZEOF_UINT_FAST8_T 0
#define H5_SIZEOF_UINT_LEAST16_T 0
#define H5_SIZEOF_UINT_LEAST32_T 0
#define H5_SIZEOF_UINT_LEAST64_T 0
#define H5_SIZEOF_UINT_LEAST8_T 0
#define H5_HAVE_DIFFTIME 1
#define H5_HAVE_FORK 1
#define H5_HAVE_GETHOSTNAME 1
#define H5_HAVE_IOCTL 1
#define H5_HAVE_LONGJMP 1
#define H5_HAVE_SIGACTION 1
#define H5_HAVE_SIGNAL 1
#define H5_HAVE_STRDUP 1
#define H5_HAVE_SYSTEM 1
#define H5_HAVE_VSNPRINTF 1
#define  H5_HAVE_IO_H 1
#define H5_HAVE_SETJMP_H 1
#define H5_HAVE_STDDEF_H 1
#define H5_HAVE_SYS_STAT_H 1
#define H5_HAVE_SYS_TIMEB 1
#define H5_HAVE_SYS_TYPES_H 1
#define H5_HAVE_WINSOCK_H 1


/* comment the following line out if the memory buffers being written to 
   disk should not be cleared before writing. */
#define H5_CLEAR_MEMORY 1

/* comment the following line out if you are not using check sum filter*/
#define H5_HAVE_FILTER_FLETCHER32 1

/* comment the following line out if you are not using shuffle filter*/
#define H5_HAVE_FILTER_SHUFFLE 1

/* comment the following two lines out if you are not using deflate(gzip) filter*/
#define H5_HAVE_FILTER_DEFLATE 1
#define H5_HAVE_ZLIB_H 1

/* comment the following two lines out if you are not using szip filter*/
/*#define H5_HAVE_SZLIB_H 1
#define H5_HAVE_FILTER_SZIP 1
*/
/* comment the following line out if your szlib library does not have an encoder */
#define H5_SZIP_CAN_ENCODE 1


#if defined(__MWERKS__) || defined(__cplusplus)
#define H5_inline inline
#else
#define H5_inline __inline
#endif

#if _MSC_VER >= 1300 /* .Net supports FUNCTION */
#define H5_HAVE_FUNCTION 1
#else
#undef H5_HAVE_FUNCTION
#endif 

#if _MSC_VER >=1400
/* visual studio 2005 doesn't support size of setvbuf to be less thn 1,This is a hacking, we would like to wait
visual studio 2005 to fix this problem.
*/

#define HDsetvbuf(F,S,M,Z) (((Z)>1)?setvbuf(F,S,M,Z):setvbuf(F,S,M,2))

#else
#define HDsetvbuf(F,S,M,Z) setvbuf(F,S,M,Z)
#endif
