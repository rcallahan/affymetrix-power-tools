////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License
// (version 2.1) as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

//
// affy/sdk/file5/File5.cpp ---
//
// $Id: File5_File.cpp,v 1.33 2009-10-20 01:37:37 rsatin Exp $
//

//
#include "file5/File5_File.h"
//
#include "file5/File5.h"
//
#include "util/AffxMultiDimensionalArray.h"
#include "util/Convert.h"
#include "util/Err.h"
#include "util/Util.h"
#include "util/Fs.h"
//
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
//

//
affx::File5_File::File5_File()
{
  m_kind_char='F';
  init();
}

affx::File5_File::~File5_File()
{
  close();
  // whine if there are objects which refer to this one.
  if (refcnt()!=0) {
#ifdef FILE5_DEBUG_PRINT
    printf("m_refcnt==%d\n",refcnt());
#endif
    //assert(0);
  }
}

affx::File5_return_t
affx::File5_File::init()
{
  File5_Group::init();
  //
  m_kind=affx::FILE5_KIND_FILE;
  m_kind_char='F';
  //
  return affx::FILE5_OK;
}

char affx::File5_File::file5_kind_char() {
  return 'F';
}

affx::File5_return_t
affx::File5_File::dump()
{
  printf("== File5_File: (%p)\n",this);
  return dump1();
}

affx::File5_return_t
affx::File5_File::dump1()
{
  refcnt_print();
  //
  printf("   m_file_name    = '%s'\n",m_file_name.c_str());
  //printf("   m_h5_file_id   = %d\n",m_h5_file_id);
  printf("   m_refcnt       = %d\n",m_refcnt);
  return affx::File5_Group::dump1();
}

//
void
affx::File5_File::setFilename(const std::string& file_name)
{
  m_file_name=file_name;
  m_name=file_name;
}

bool
affx::File5_File::isHdf5file(const std::string& file_name)
{
  std::string tmp_unc_path=file_name;
  Fs::convertToUncPathInPlace(tmp_unc_path,10);

  if (Fs::isReadable(tmp_unc_path) &&
      H5Fis_hdf5(tmp_unc_path.c_str())==1) {
    return true;
  }
  //
  return false;
}


//
affx::File5_return_t
affx::File5_File::open(const std::string& file_name,int flags)
{
  //int rv;
  //struct stat stat_buf;
  std::string tmp_unc_path=file_name;
  Fs::convertToUncPathInPlace(tmp_unc_path, 10);

  if (refcnt()!=0) {
    FILE5_ABORT("file still has "+ToStr(refcnt())+" open references.");
  }
  
  setFilename(file_name);
  m_flags=flags;
  m_readonly=true;

  // replace means get rid of the orginal and implies create
  if ((flags&FILE5_REPLACE)==FILE5_REPLACE) {
#ifdef FILE5_DEBUG_PRINT
    printf("### file5: file replace: '%s'\n",m_file_name.c_str());
#endif
    m_h5_obj=H5Fcreate(tmp_unc_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if ( m_h5_obj < 1 ) {
      Verbose::out(1, "H5Fcreate failed on absolute path, trying relative..." + file_name);
      m_h5_obj=H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if ( m_h5_obj >= 0 ) {
        Verbose::out(1, "H5Fcreate ok for relative path: " + file_name);
      }
    }    
    FILE5_CHECKID(m_h5_obj,"could not replace "+FS_QUOTE_PATH(file_name));
    m_readonly=false;
  }
  //
  else if (((flags&FILE5_OPEN)==FILE5_OPEN)&&(isHdf5file(m_file_name)==1)) {
    // printf("### file5: file open: '%s'\n",m_file_name.c_str());
    // check to see if it exists
    if (!Fs::isReadable(tmp_unc_path)) {
      FILE5_ABORT("File does not exist or is not readable: "+FS_QUOTE_PATH(tmp_unc_path));
    }
    //
    int h5_open_flags;
    if ((flags&affx::FILE5_RO)==affx::FILE5_RO) {
      h5_open_flags=H5F_ACC_RDONLY;
      m_readonly=true;
    }
    else {
      h5_open_flags=H5F_ACC_RDWR;
      m_readonly=false;
    }
    m_h5_obj=H5Fopen(tmp_unc_path.c_str(),h5_open_flags,H5P_DEFAULT);
    FILE5_CHECKID(m_h5_obj,"could not open: "+FS_QUOTE_PATH(tmp_unc_path));
  }
  else if ((flags&FILE5_CREATE)==FILE5_CREATE) {
#ifdef FILE5_DEBUG_PRINT
    printf("### file5: file create: '%s'\n",tmp_unc_path.c_str());
#endif
    m_h5_obj=H5Fcreate(tmp_unc_path.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    FILE5_CHECKID(m_h5_obj,"could not create: "+FS_QUOTE_PATH(tmp_unc_path));
    //
    m_readonly=false;
  }
  else if ((flags&FILE5_OPEN)==FILE5_OPEN) {
    printf("### file5: file open: unable to open: '%s'\n",tmp_unc_path.c_str());
    return affx::FILE5_ERR;
  }
  else {
    FILE5_ABORT("Bad set of flags passed to File5_File::open()");
  }

  // we dont bump the refcnt, we loop back the parent to ourselves
  setParent(this);
  setFile(this);
  m_state=affx::FILE5_STATE_OPEN;

  //
#ifdef FILE5_DEBUG_PRINT
  printf("### File5_File::open()==%d\n",m_h5_obj);
#endif

  return affx::FILE5_OK;
}

bool
affx::File5_File::is_open()
{
  // using "0" and "1" as this returns a boolean.
  if (m_h5_obj==-1) {
    return false;
  }
  else {
    return true;
  }
}

affx::File5_return_t
affx::File5_File::close()
{
  flush();

  // warn about unclosed references.
  if (refcnt()!=0) {
    printf("File5_File::close('%s'): refcnt==%d\n",m_file_name.c_str(),refcnt());
    // abort before deallocating resources for debugging.
    // FILE5_ABORT("File5_File::close(): refcnt!=0");
    // return 0;
  }

  // we will close the file ourselves, with H5Fclose.
  // let the group have its turn.
  hid_t tmp_h5_obj=m_h5_obj;
  m_h5_obj=-1;
  affx::File5_Group::close();
  //
  if (tmp_h5_obj!=-1) {
    m_h5_status=H5Fclose(tmp_h5_obj);
    FILE5_CHECKRV(m_h5_status,"H5Fclose failed.");
  }
  //
  return affx::FILE5_OK;
}

affx::File5_return_t
affx::File5_File::flush()
{
  // for debugging to check that flush is called.
  // printf("FLUSH ==========\n"); sleep(2);
  if (m_h5_obj!=-1) {
    m_h5_status=H5Fflush(m_h5_obj,H5F_SCOPE_GLOBAL);
    FILE5_CHECKRV(m_h5_status,"H5Fflush failed.");
  }
  return affx::FILE5_OK;
}

// I think this is a better signature. - jhg
// bool affx::File5_File::equivalent(
// const std::string& strFileName1,
// const std::string& strTsvName1,
// const std::string& strFileName2,
// const std::string& strTsvName2,
// std::set<std::string>& setIgnore, 
// double dEpsilon, 
// double dCorrelationCutoff)
bool affx::File5_File::equivalent(	const std::string& strFileName1, 
                                        const std::string& strFileName2, 
                                        const std::string& strGroupName, 
                                        const std::string& strTsvName, std::set<std::string>& setIgnore, 
                                        double dEpsilon, 
                                        double dCorrelationCutoff,
                                        bool bAllowNegation,
                                        bool flagNaNNumDiff)
{
	Verbose::out(1, strFileName1);
	Verbose::out(1, strFileName2);
	Verbose::out(1, "Comparing " + strGroupName + "." + strTsvName);
	bool bSuccessful = true;
	double d1 = 0;
	double d2 = 0;
	float f1 = 0;
	float f2 = 0;
	int i1 = 0;
	int i2 = 0;
	char c1 = 0;
	char c2 = 0;
	std::string str1;
	std::string str2;

	if (!isHdf5file(strFileName1)) {Verbose::out(1, "File: " + strFileName1 + " is not an HDF5 file."); return false;}
	if (!isHdf5file(strFileName2)) {Verbose::out(1, "File: " + strFileName2 + " is not an HDF5 file."); return false;}

	int iCount1 = affx::File5_Tsv::getFileTsvLineCount(strFileName1, strGroupName, strTsvName);
	int iCount2 = affx::File5_Tsv::getFileTsvLineCount(strFileName2, strGroupName, strTsvName);
	if (iCount1 != iCount2) {Verbose::out(1, "Tsv Line Counts do not match for Group: " + strGroupName + ", Tsv: " + strTsvName); return false;}
	if (iCount1<0 && iCount2<0) {Verbose::out(1, "Missing from both files, skipping test for Group: " + strGroupName + ", Tsv: " + strTsvName); return true;}  // files match if both files do not have the Group and/or TSV Column

	affx::File5_File file1;
	affx::File5_File file2;
	// affx::File5_File* file1=new File5_File();
	affx::File5_Group* group1 = NULL;
	affx::File5_Group* group2 = NULL;
	affx::File5_Tsv* tsv1 = NULL;
	affx::File5_Tsv* tsv2 = NULL;
	file1.open(strFileName1, affx::FILE5_OPEN_RO);
	file2.open(strFileName2, affx::FILE5_OPEN_RO);
	group1 = file1.openGroup(strGroupName, affx::FILE5_OPEN_RO);
	group2 = file2.openGroup(strGroupName, affx::FILE5_OPEN_RO);
	tsv1 = group1->openTsv(strTsvName, affx::FILE5_OPEN_RO);
	tsv2 = group2->openTsv(strTsvName, affx::FILE5_OPEN_RO);
	/// @todo the above is too complicated - jhg
	//  They should both be this:
	// tsv1=file1->openTsv(strTsvName1,affx::FILE5_OPEN_RO)
	int iColumnCount1 = tsv1->getColumnCount(0);
	int iColumnCount2 = tsv2->getColumnCount(0);
	if (iColumnCount1 != iColumnCount2) {Verbose::out(1, "Tsv Column Counts do not match for Group: " + strGroupName + ", Tsv: " + strTsvName); bSuccessful = false;}
	else
	{
		for (int iColumnIndex = 0; (iColumnIndex < iColumnCount1); iColumnIndex++)
		{
			tsv1->close();
			tsv2->close();
			delete tsv1;
			delete tsv2;
			tsv1 = group1->openTsv(strTsvName, affx::FILE5_OPEN_RO);
			tsv2 = group2->openTsv(strTsvName, affx::FILE5_OPEN_RO);
			double dMaxDifference = 0;
			AffxMultiDimensionalArray<double> mx(iCount1, 2);
			std::string strColumnName1;
			std::string strColumnName2;
			tsv1->getColumnName(0, iColumnIndex, &strColumnName1);
			tsv2->getColumnName(0, iColumnIndex, &strColumnName2);
			if (strColumnName1 != strColumnName2) {Verbose::out(1, "Tsv Column Names do not match for Group: " + strGroupName + ", Tsv: " + strTsvName + ", Column: " + ToStr(iColumnIndex)); bSuccessful = false; continue;}
					
			std::string strColumnName = strGroupName + "." + strTsvName + "." + strColumnName1;
			if (setIgnore.find(strColumnName) != setIgnore.end()) {
				Verbose::out(1, "Skipping test for column : " + strColumnName);
				continue;
			}

			affx::File5_dtype_t iColumnType1 = tsv1->getColumnDtype(0, iColumnIndex);
			affx::File5_dtype_t iColumnType2 = tsv2->getColumnDtype(0, iColumnIndex);
			if (iColumnType1 != iColumnType2) {Verbose::out(1, "Tsv Column Types do not match for Group: " + strGroupName + ", Tsv: " + strTsvName + ", Column: " + strColumnName1); bSuccessful = false; continue;}

			bool bNaNDiff = false;
			int iRowIndex = 0;
			while ((tsv1->nextLine() == affx::FILE5_OK) && (tsv2->nextLine() == affx::FILE5_OK))
			{
				if (iColumnType1 == affx::FILE5_DTYPE_STRING)
				{
					tsv1->get(0, iColumnIndex, &str1);
					tsv2->get(0, iColumnIndex, &str2);
					if (str1 != str2) {Verbose::out(1, "Tsv Column Values do not match for Group: " + strGroupName + ", Tsv: " + strTsvName + ", Column: " + strColumnName1 + ", RowIndex = " + ToStr(iRowIndex) + ", Value1: " + str1 + ", Value2: " + str2); bSuccessful = false;}
				}
				else
				{
					switch (iColumnType1)
					{
					case affx::FILE5_DTYPE_DOUBLE:
						tsv1->get(0, iColumnIndex, &d1);
						tsv2->get(0, iColumnIndex, &d2);
						break;
					case affx::FILE5_DTYPE_FLOAT:
						tsv1->get(0, iColumnIndex, &f1);
						tsv2->get(0, iColumnIndex, &f2);
						d1 = f1;
						d2 = f2;
						break;
					case affx::FILE5_DTYPE_INT:
						tsv1->get(0, iColumnIndex, &i1);
						tsv2->get(0, iColumnIndex, &i2);
						d1 = i1;
						d2 = i2;
						break;
					case affx::FILE5_DTYPE_CHAR:
						tsv1->get(0, iColumnIndex, &c1);
						tsv2->get(0, iColumnIndex, &c2);
						d1 = c1;
						d2 = c2;
						break;
					default:
						APT_ERR_ABORT("Unhandled case.")
						break;
					}
					if(bAllowNegation)
					{
						mx.set(iRowIndex, 0, fabs(d1));
						mx.set(iRowIndex, 1, fabs(d2));
						// Max can yield invalid results if either operand is NaN
						if ((d1 != d1) || (d2 != d2))
						{
							// If flag NaN-num diff is true and not both d1 and d2 are NaN but one is, set flag.
							if (flagNaNNumDiff && !((d1 != d1) && (d2 != d2))) { bNaNDiff = true; }
						}
						else
						{
							dMaxDifference = Max(dMaxDifference, fabs(fabs(d1) - fabs(d2)));
						}
					}
					else
					{
						mx.set(iRowIndex, 0, d1);
						mx.set(iRowIndex, 1, d2);
						// Max can yield invalid results if either operand is NaN
						if (((d1 != d1) || (d2 != d2)))
						{
							// If flag NaN-num diff is true and not both d1 and d2 are NaN but one is, set flag.
							if (flagNaNNumDiff && !((d1 != d1) && (d2 != d2))) {bNaNDiff = true; }
						}
						else
						{
							dMaxDifference = Max(dMaxDifference, fabs(d1 - d2));
						}
					}
				}
				iRowIndex++;
			}
			double dCorrelation = mx.corr();					
//			Verbose::out(1, strGroupName + "." + strTsvName + "." + strColumnName1 + ", MaxDifference = " + ToStr(dMaxDifference) + ", Correlation:: " + ToStr(dCorrelation)); 
			if (dMaxDifference > dEpsilon || bNaNDiff)
			{
				if ((fabs(dCorrelation) < dCorrelationCutoff) || (dCorrelationCutoff == 1) || bNaNDiff)
				{
#ifdef WIN32
					if (bNaNDiff) { dMaxDifference = ::sqrt(-1.); } // Sets to -1.#IND000000000000 and not 1.#QNAN like numeric_limits=<double>::quiet_NaN()
#else
					if (bNaNDiff) { dMaxDifference = numeric_limits<double>::quiet_NaN(); }
#endif
					Verbose::out(1, "Tsv Column Values are out of spec. for Group: " + strGroupName + ", Tsv: " + strTsvName + ", Column: " + strColumnName1 + ", MaxDifference = " + ToStr(dMaxDifference) + ", Correlation: " + ::getDouble(dCorrelation, 12)); 
					bSuccessful = false;
				}
			}
		}	
	}
	tsv1->close();
	tsv2->close();
	delete tsv1;
	delete tsv2;

	group1->close();
	group2->close();
	delete group1;
	delete group2;

	file1.close();
	file2.close();

	if (bSuccessful) {Verbose::out(1, "TSVs are equivalent.");}

	return bSuccessful;
}

