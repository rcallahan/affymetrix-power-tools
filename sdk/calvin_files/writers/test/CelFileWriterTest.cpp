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

#include "calvin_files/writers/test/CelFileWriterTest.h"
//
#include "calvin_files/writers/src/CalvinCelFileWriter.h"
//

using namespace std;
using namespace affymetrix_calvin_io;

CPPUNIT_TEST_SUITE_REGISTRATION( CelFileWriterTest );

void CelFileWriterTest::setUp() {}

void CelFileWriterTest::tearDown() {}

void CelFileWriterTest::testCreation()
{
	CelFileData fHdr("cel_file");
	CelFileWriter* w = new CelFileWriter(fHdr);
	CPPUNIT_ASSERT(1);
	delete w;
}

void CelFileWriterTest::WriteTest()
{
	std::string celFileName = "cel_file";
	CelFileData writerData(celFileName);
	writerData.SetIntensityCount(10);
	writerData.SetStdDevCount(10);
	writerData.SetPixelCount(10);
	writerData.SetOutlierCount(10);
	writerData.SetMaskCount(10);

    writerData.SetArrayType(L"arraytype");
    writerData.SetMasterFileName(L"masterfile");
    writerData.SetLibraryPackageName(L"libpackage");

	CelFileWriter writer(writerData);

	FloatVector intensities;
	intensities.push_back(10.0);
	intensities.push_back(11.0);
	intensities.push_back(12.0);
	intensities.push_back(13.0);
	writer.WriteIntensities(intensities);
	CPPUNIT_ASSERT(1);

	FloatVector stdDevs;
	stdDevs.push_back(10.0);
	stdDevs.push_back(11.0);
	stdDevs.push_back(12.0);
	stdDevs.push_back(13.0);
	writer.WriteStdDevs(stdDevs);
	CPPUNIT_ASSERT(1);

	Int16Vector pixels;
	pixels.push_back(10);
	pixels.push_back(11);
	pixels.push_back(12);
	pixels.push_back(13);
	writer.WritePixels(pixels);
	CPPUNIT_ASSERT(1);

	CelFileData readerData;
	CelFileReader reader;
	reader.SetFilename(celFileName);
	CPPUNIT_ASSERT_NO_THROW(reader.Read(readerData));

	// Read some data
	CPPUNIT_ASSERT(readerData.GetFileHeader()->GetFilename() == celFileName);
	CPPUNIT_ASSERT(readerData.HasNumPixels());
	CPPUNIT_ASSERT(readerData.HasStdev());

    
    CPPUNIT_ASSERT(readerData.GetArrayType() == L"arraytype");
    CPPUNIT_ASSERT(readerData.GetMasterFileName() == L"masterfile");
    CPPUNIT_ASSERT(readerData.GetLibraryPackageName() == L"libpackage");
}
