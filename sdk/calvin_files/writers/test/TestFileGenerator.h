////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#ifndef __TESTFILEGENERATOR_H_
#define __TESTFILEGENERATOR_H_

/*! \file TestFileGenerator.h This defines a class that generates test data files. 
 *  This is not really a test fixture.  Once the files are
 *	generated they still have to be copied to data/data or
 *	parsers/data or fusion/data folders.  To generate a
 *	file, the method that generates it must be uncommented
 *	in GenerateTestFiles.
 */

#include "calvin_files/data/src/CELData.h"
#include "calvin_files/data/src/DATData.h"
#include "calvin_files/data/src/GenericDataHeader.h"
//
#include <cppunit/extensions/HelperMacros.h>
//

class TestFileGenerator : public CPPUNIT_NS::TestFixture  
{
	CPPUNIT_TEST_SUITE( TestFileGenerator );

	CPPUNIT_TEST ( GenerateTestFiles );

	CPPUNIT_TEST_SUITE_END();

public:

	void GenerateTestFiles();

private:

// METHODS TO GENERATE FILES.
	void WriteOutGenericDATDataFile1UsingRawWrites();
	void WriteOutGenericDATDataFile1UsingGenericWriter();
	void WriteOutGenericDATDataFileNoGrid();
	void WriteOutGenericDATDataFileWithGrid();
	void WriteOutGenericDataFileWithAllColumnTypes();
	void WriteSmallCelFileNoOutlierNoMask();
	void WriteSmallCelFileNoStdev();
	void WriteSmallCelFile();
	//void WriteSmallUint16CelFile();	// Requires changes in CelFileWriter and CELData.
	void WriteLargeCelFile();
	void WriteCelFileWithADataSetWithZeroRows();
	void WriteSmallDatFile();
	void WriteLargeDatFile();
	void WriteActualSizeExpressionCDFFile();
	void WriteSmallExpressionCDFFile();
	void WriteSmallQCCDFFile();
	void WriteSmallDatFileWithGridAndSubgrids();
	void WriteSmallCelFileWithAPartialDatHeaderTest();
	void WriteSmallCelFileWithAFullDatHeaderTest();
	void WriteSmallDatFileWithGridAndSubgridsAndParameters();
	void WriteSmallDatFileWithReservedStringParameters();

// HELPER METHODS
	void AddStandardGenericDataHeader(affymetrix_calvin_io::GenericDataHeader& gdh);
	void WriteRemaingSmallCelFileWithGridParameters(affymetrix_calvin_io::CelFileData& data);
	void WriteDatFile(const std::string& name, const std::wstring& type, int32_t rows, int32_t cols, bool showProgress);
	void AddGridAndSubgrids(affymetrix_calvin_io::DATData& data, float increment, int32_t subgridCnt);
	void WriteExpressionCDFFile(const std::string& filename, u_int32_t probeSetCnt);
	void WriteQCCDFFile(const std::string& filename, u_int32_t probeSetCnt);

};

#endif // __TESTFILEGENERATOR_H_
