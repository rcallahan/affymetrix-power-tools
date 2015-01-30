////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////
//
// CalvinLite/CL_util.h ---
//
// $Id: CL_util.h,v 1.1 2009-10-27 16:53:53 harley Exp $
//

//
#include "calvinlite/CalvinLite.h"
//
#include <stdexcept>
#include <string>

// What CalvinLite Objects will throw when it has an error, unless "throwOnErr(false)" is called.
class CalvinLiteException : public std::runtime_error {
public:
  /// the numeric error thrown.
  CL_Err_t m_err_num;
  /// the text of the error message.
  std::string m_err_msg;
  /// The path of the file which cause the error. (so it can be reported.)
  std::string m_err_path;

  CalvinLiteException(CL_Err_t err_num,
                      const std::string& err_msg,
                      const std::string& err_path
                      ) : std::runtime_error("CalvinLite")
  {
    m_err_num=err_num;
    m_err_msg=err_msg;
    m_err_path=err_path;
  };
  // @todo: why is the "throw ()" required?
  ~CalvinLiteException() throw () { };
};


// see the note in CL_gen_uuid about setting the seed.
extern int CL_gen_uuid_seed;

//
std::string CL_gen_uuid();

//
std::string CL_gen_timestr();
std::string CL_gen_timestr(time_t when);

//
extern const char* CL_error_strings[];
std::string CL_get_err_string(CL_Err_t err);

//
CL_TypeCode_t CL_mimeStrToCode(const std::string& mimeStr);
std::string CL_codeToMimeStr(CL_TypeCode_t tcode);
int CL_codeToBytelen(CL_TypeCode_t tcode);

//
int   CL_mem_read_u1B(void* ptr);
int   CL_mem_read_i1B(void* ptr);
int   CL_mem_read_u2B(void* ptr);
int   CL_mem_read_i2B(void* ptr);
int   CL_mem_read_u4B(void* ptr);
int   CL_mem_read_i4B(void* ptr);
float CL_mem_read_float(void* ptr);

//
void CL_mem_write_i1B(void* ptr,int val);
void CL_mem_write_i2B(void* ptr,int val);
void CL_mem_write_i4B(void* ptr,int val);
void CL_mem_write_float(void* ptr,float val);

//
std::string CL_basename(const std::string& path);

